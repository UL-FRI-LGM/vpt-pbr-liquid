// #package js/main

// #include utils
// #include EventEmitter.js
// #include WebGL.js
// #include Draggable.js

// #include ../html/TransferFunctionWidget.html
// #include ../html/TransferFunctionWidgetBumpHandle.html
// #include ../css/TransferFunctionWidget.css

class TransferFunctionWidget extends EventEmitter {

constructor(options) {
    super();
    this._removeHandleByKey = this._removeHandleByKey.bind(this);

    this._onColorChange = this._onColorChange.bind(this);

    Object.assign(this, {
        _width                  : 256,
        _height                 : 256,
        _transferFunctionWidth  : 256,
        _transferFunctionHeight : 256,
        scaleSpeed              : 0.003
    }, options);

    this._$html = DOMUtils.instantiate(TEMPLATES.TransferFunctionWidget);
    this._$colorPicker          = this._$html.querySelector('[name="color"]');
    this._$alphaPicker          = this._$html.querySelector('[name="alpha"]');
    this._$channelContribPicker = this._$html.querySelector('[name="channel-contribution"]');
    this._$addBumpButton        = this._$html.querySelector('[name="add-bump"]');
    this._$removeSelectedBump   = this._$html.querySelector('[name=remove-selected-bump]')
    this._$removeAllBumps       = this._$html.querySelector('[name=remove-all-bumps]')
    this._$loadButton           = this._$html.querySelector('[name="load"]');
    this._$saveButton           = this._$html.querySelector('[name="save"]');

    this._canvas = this._$html.querySelector('canvas');
    this._canvas.width = this._transferFunctionWidth;
    this._canvas.height = this._transferFunctionHeight;
    this.resize(this._width, this._height);

    this._gl = this._canvas.getContext('webgl2', {
        depth                 : false,
        stencil               : false,
        antialias             : false,
        preserveDrawingBuffer : true
    });
    const gl = this._gl;

    gl.enable(gl.BLEND);
    gl.blendFunc(gl.SRC_ALPHA, gl.ONE_MINUS_SRC_ALPHA);

    this._clipQuad = WebGL.createClipQuad(gl);
    this._program = WebGL.buildPrograms(gl, {
        drawTransferFunction: SHADERS.drawTransferFunction
    }, MIXINS).drawTransferFunction;
    const program = this._program;
    gl.useProgram(program.program);
    gl.bindBuffer(gl.ARRAY_BUFFER, this._clipQuad);
    gl.enableVertexAttribArray(program.attributes.aPosition);
    gl.vertexAttribPointer(program.attributes.aPosition, 2, gl.FLOAT, false, 0, 0);

    this._bumps = [];
    this._$addBumpButton.addEventListener('click', () => {
        this.addBump();
    });
    this._$removeSelectedBump.addEventListener('click', () => {
        this._removeSelectedBump();
    });
    this._$removeAllBumps.addEventListener('click', () => {
        this._removeAllBumps();
    });

    this._$colorPicker.addEventListener('input', this._onColorChange);
    this._$alphaPicker.addEventListener('input', this._onColorChange);

    this._channelContribution = 1.0;
    this._$channelContribPicker.addEventListener('input', this._onColorChange);

    this._$loadButton.addEventListener('click', () => {
        CommonUtils.readTextFile(data => {
            this.loadFromJson(JSON.parse(data));
        });
    });

    this._$saveButton.addEventListener('click', () => {
        CommonUtils.downloadJSON(this._bumps, 'TransferFunction.json');
    });
}

destroy() {
    const gl = this._gl;
    gl.deleteBuffer(this._clipQuad);
    gl.deleteProgram(this._program.program);
    DOMUtils.remove(this._$html);
}

resize(width, height) {
    this._canvas.style.width = width + 'px';
    this._canvas.style.height = height + 'px';
    this._width = width;
    this._height = height;
}

resizeTransferFunction(width, height) {
    this._canvas.width = width;
    this._canvas.height = height;
    this._transferFunctionWidth = width;
    this._transferFunctionHeight = height;
    const gl = this._gl;
    gl.viewport(0, 0, gl.drawingBufferWidth, gl.drawingBufferHeight);
}

render() {
    const gl = this._gl;
    const program = this._program;

    gl.clear(gl.COLOR_BUFFER_BIT);
    this._bumps.forEach(bump => {
        gl.uniform2f(program.uniforms['uPosition'], bump.position.x, bump.position.y);
        gl.uniform2f(program.uniforms['uSize'], bump.size.x, bump.size.y);
        gl.uniform4f(program.uniforms['uColor'], bump.color.r, bump.color.g, bump.color.b, bump.color.a);
        gl.drawArrays(gl.TRIANGLE_FAN, 0, 4);
    });
}

loadFromJson(json) {
    this._bumps = json;
    this.render();
    this._rebuildHandles();
    this.trigger('change');
}

addBump(options) {
    const bumpIndex = this._bumps.length;
    const newBump = {
        position: {
            x: 0.5,
            y: 0.5
        },
        size: {
            x: 0.2,
            y: 0.2
        },
        color: {
            r: 1,
            g: 0,
            b: 0,
            a: 1
        }
    };
    this._bumps.push(newBump);
    this._addHandle(bumpIndex);
    this.selectBump(bumpIndex);
    this.render();
    this.trigger('change');
}

_addHandle(index) {
    const $handle = DOMUtils.instantiate(TEMPLATES.TransferFunctionWidgetBumpHandle);
    this._$html.querySelector('.widget').appendChild($handle);
    DOMUtils.data($handle, 'index', index);

    const left = this._bumps[index].position.x * this._width;
    const top = (1 - this._bumps[index].position.y) * this._height;
    $handle.style.left = Math.round(left) + 'px';
    $handle.style.top = Math.round(top) + 'px';

    new Draggable($handle, $handle.querySelector('.bump-handle'));
    $handle.addEventListener('draggable', e => {
        const x = Math.max(0, Math.min(e.currentTarget.offsetLeft / this._width, 1));
        const y = Math.min(1, Math.max(0, 1 - (e.currentTarget.offsetTop / this._height)));
        const i = parseInt(DOMUtils.data(e.currentTarget, 'index'));
        this._bumps[i].position.x = x;
        this._bumps[i].position.y = y;
        this.render();
        this.trigger('change');
    });
    $handle.addEventListener('mousedown', e => {
        const i = parseInt(DOMUtils.data(e.currentTarget, 'index'));
        this.selectBump(i);
    });
    $handle.addEventListener('mousewheel', e => {
        const amount = e.deltaY * this.scaleSpeed;
        const scale = Math.exp(-amount);
        const i = parseInt(DOMUtils.data(e.currentTarget, 'index'));
        this.selectBump(i);
        if (e.shiftKey) {
            this._bumps[i].size.y *= scale;
        } else {
            this._bumps[i].size.x *= scale;
        }
        this.render();
        this.trigger('change');
    });
    $handle.addEventListener('mouseover', () => {
        // console.log("hello");
        document.addEventListener('keydown', this._removeHandleByKey);
    });
    $handle.addEventListener('mouseout', () => {
        // console.log("bye");
        document.removeEventListener('keydown', this._removeHandleByKey);
    });
}

_removeHandleByKey(e) {
    if (e.key == 'Delete')
        this._removeSelectedBump();
}

_removeSelectedBump() {
    this._removeHandle(this.getSelectedBumpIndex());
}

_removeAllBumps() {
    this._bumps = [];
    this._rebuildHandles();
    this.render();
    this.trigger('change');
}

_removeHandle(index) {
    const handles = this._$html.querySelectorAll('.bump');
    handles.forEach(handle => {
        const i = parseInt(DOMUtils.data(handle, 'index'));
        if (i == index) {
            this._bumps.splice(i, 1);
        }
    });
    this._rebuildHandles();
    this.render();
    this.trigger('change');
}

_rebuildHandles() {
    const handles = this._$html.querySelectorAll('.bump');
    handles.forEach(handle => {
        DOMUtils.remove(handle);
    });
    for (let i = 0; i < this._bumps.length; i++) {
        this._addHandle(i);
    }
}

selectBump(index) {
    const handles = this._$html.querySelectorAll('.bump');
    handles.forEach(handle => {
        const i = parseInt(DOMUtils.data(handle, 'index'));
        if (i === index) {
            handle.classList.add('selected');
        } else {
            handle.classList.remove('selected');
        }
    });

    const color = this._bumps[index].color;
    this._$colorPicker.value = CommonUtils.rgb2hex(color);
    this._$alphaPicker.value = color.a;
}

getSelectedBumpIndex() {
    const handles = this._$html.querySelectorAll('.bump');
    let idx = -1;
    handles.forEach(handle => {
        let i = parseInt(DOMUtils.data(handle, 'index'));
        if (handle.classList.contains('selected')) {
            idx = i;
        }
    });
    return idx;
}

getTransferFunction() {
    return this._canvas;
}

_onColorChange() {
    const $selectedBump = this._$html.querySelector('.bump.selected');
    if ($selectedBump !== null) {
        const i = parseInt(DOMUtils.data($selectedBump, 'index'));
        const color = CommonUtils.hex2rgb(this._$colorPicker.value);
        const alpha = parseFloat(this._$alphaPicker.value);
        this._bumps[i].color.r = color.r;
        this._bumps[i].color.g = color.g;
        this._bumps[i].color.b = color.b;
        this._bumps[i].color.a = alpha;
    }
    this._channelContribution = parseFloat(this._$channelContribPicker.value);

    this.render();
    this.trigger('change');
}

appendTo(object) {
    object.appendChild(this._$html);
}

disableButtons(value) {
    this._$addBumpButton.disabled = value;
    this._$removeSelectedBump.disabled = value;
    this._$removeAllBumps.disabled = value;
    this._$loadButton.disabled = value;
    this._$saveButton.disabled = value;
}

}
