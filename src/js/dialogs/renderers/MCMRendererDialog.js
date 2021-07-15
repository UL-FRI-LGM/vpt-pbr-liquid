// #package js/main

// #include ../AbstractDialog.js
// #include ../../TransferFunctionWidget.js

// #include ../../../uispecs/renderers/MCMRendererDialog.json

// #include ../../ui

class MCMRendererDialog extends AbstractDialog {

constructor(renderer, options) {
    super(UISPECS.MCMRendererDialog, options);

    this._renderer = renderer;

    this._handleChange = this._handleChange.bind(this);
    this._handleTFChange = this._handleTFChange.bind(this);

    this._binds.extinction.addEventListener('input', this._handleChange);
    this._binds.albedo.addEventListener('change', this._handleChange);
    this._binds.bias.addEventListener('change', this._handleChange);
    this._binds.ratio.addEventListener('change', this._handleChange);
    this._binds.bounces.addEventListener('input', this._handleChange);
    this._binds.steps.addEventListener('input', this._handleChange);
    this._binds.minCutPlanes.addEventListener("input", this._handleChange);
    this._binds.maxCutPlanes.addEventListener("input", this._handleChange);
    this._binds.viewCutDistance.addEventListener('input', this._handleChange);
    this._binds.maxContribution.addEventListener('change', this._handleChange);
    this._binds.origData.addEventListener('change', this._handleChange);
    this._binds.origVsSeg.addEventListener('change', this._handleChange);
    this._binds.bilateral.addEventListener('change', this._handleChange);
    this._binds.bilateralGradient.addEventListener('change', this._handleChange);
    this._binds.bilateralSigma.addEventListener('input', this._handleChange);
    this._binds.bilateralBSigma.addEventListener('input', this._handleChange);
    this._binds.bilateralMSize.addEventListener('input', this._handleChange);

    this._tfwidgets = [];
    this._customTf = [];
    for (let i = 0; i < 4; i++) {
        this._tfwidgets[i] = new TransferFunctionWidget();
        let panel = new Panel();
        panel.add(this._tfwidgets[i]);
        this._customTf[i] = this._createUIElementsForCustomTf(i);
        panel.add(this._customTf[i]);
        this._binds.tftabs.add(""+i, panel);
        this._tfwidgets[i].addEventListener('change', () => {this._handleTFChange(i)});
        this._customTf[i].addEventListener('change', () => {this._handleCustomTfChange(i)});
    }
}

destroy() {
    if (this._tfwidget)
        this._tfwidget.destroy();
    super.destroy();
}

_updateTFWidgets(number) {
    for (let i = 1; i <= 4; i++) {
        if (i <= number)
            this._binds.tftabs._binds.headers.children[i].hidden = false;
        else
            this._binds.tftabs._binds.headers.children[i].hidden = true;
    }
}

_handleChange() {
    const extinction = this._binds.extinction.getValue();
    const albedo     = this._binds.albedo.getValue();
    const bias       = this._binds.bias.getValue();
    const ratio      = this._binds.ratio.getValue();
    const bounces    = this._binds.bounces.getValue();
    const steps      = this._binds.steps.getValue();
    const minCutPlaneValues = this._binds.minCutPlanes.getValue();
    const maxCutPlaneValues = this._binds.maxCutPlanes.getValue();
    const viewCutDistance   = this._binds.viewCutDistance.getValue();
    const maxContrib = this._binds.maxContribution.checked;
    const origData   = this._binds.origData.checked;
    const origVsSeg  = this._binds.origVsSeg.getValue();
    const bilateral  = this._binds.bilateral.checked;
    const bilateralGradient  = this._binds.bilateralGradient.checked;
    const bilateralSigma      = this._binds.bilateralSigma.getValue();
    const bilateralBSigma     = this._binds.bilateralBSigma.getValue();
    const bilateralMSize      = this._binds.bilateralMSize.getValue();

    this._renderer.absorptionCoefficient = extinction * (1 - albedo);
    this._renderer.scatteringCoefficient = extinction * albedo;
    this._renderer.scatteringBias = bias;
    this._renderer.majorant = extinction * ratio;
    this._renderer.maxBounces = bounces;
    this._renderer.steps = steps;
    this._renderer.minCutPlaneValues = minCutPlaneValues;
    this._renderer.maxCutPlaneValues = maxCutPlaneValues;
    this._renderer.viewCutDistance = viewCutDistance;
    this._renderer.maxContribution = maxContrib;
    this._renderer.origData = origData;
    this._renderer.origVsSeg = origVsSeg;
    this._renderer.bilateral = bilateral;
    this._renderer.bilateralGradient = bilateralGradient;
    this._renderer.bilateralSigma = bilateralSigma;
    this._renderer.bilateralBSigma = bilateralBSigma;
    this._renderer.bilateralMSize = bilateralMSize;

    this._renderer.reset();
}

_handleTFChange(id) {
    this._renderer.setTransferFunction(this._tfwidgets[id].getTransferFunction(), id);

    this._renderer._channelContributions.x = this._tfwidgets[0]._channelContribution;
    this._renderer._channelContributions.y = this._tfwidgets[1]._channelContribution;
    this._renderer._channelContributions.z = this._tfwidgets[2]._channelContribution;
    this._renderer._channelContributions.w = this._tfwidgets[3]._channelContribution;

    this._renderer.reset();
}

_createUIElementsForCustomTf(id) {
    var field = new Field({
        label: "Custom Transfer Function:"
    });
    var checkbox = new Checkbox({
        checked: false
    });
    checkbox.addEventListener('change', () => {this._handleCustomTfChange(id)});
    field.add(checkbox);
    return field;
}

_handleCustomTfChange(id) {
    var checkbox = this._customTf[id]._content;
    if (checkbox.isChecked()) {
        console.log("custom transfer function");
        this._renderer._customTf[id] = true;
        this._tfwidgets[id].disableButtons(true);
        // todo: set to custom transfer function here
        
        //this._tfwidgets[id].createCustomTransferFunction();
        //this._tfwidgets[id]._removeAllBumps();
        //this._renderer.setTransferFunction(this._tfwidgets[id].getTransferFunction(), id);
        //this.renderer.setTransferFunction();
    } else {
        console.log("normal transfer function");
        this._renderer._customTf[id] = false;
        this._tfwidgets[id].disableButtons(false);
        //this._renderer.setTransferFunction(this._tfwidgets[id].getTransferFunction(), id);
    }
    let tfArrays = this._renderer._changeCustomTf();
    tfArrays.forEach((tfArray, index) => {
        if (tfArray.length === 0) {
            // empty array, meaning no volume is selected
            return;
        }

        console.log("in foreach");
        const imgData = new ImageData(Uint8ClampedArray.from(tfArray), 256, 256);
        const canv = document.createElement('canvas');
        canv.width = 256;
        canv.height = 256;
        const ctx = canv.getContext('2d');
        ctx.putImageData(imgData, 0, 0);
        let imgDataUrl = canv.toDataURL();

        this._renderer.setTransferFunction(canv, index);

        console.log(imgDataUrl);
        this._tfwidgets[index]._canvas.style.backgroundImage = "none";
        this._tfwidgets[index]._canvas.style.backgroundImage = 'url('+imgDataUrl+')';
    });
}

_handleScaleChange(dimensions) {
    this._binds.scale.setValue(new Vector(dimensions.x, dimensions.y, dimensions.z));
}

}
