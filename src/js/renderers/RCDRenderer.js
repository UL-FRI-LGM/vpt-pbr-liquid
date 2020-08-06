// #package js/main

// #include ../WebGL.js
// #include AbstractRenderer.js
// #include ../LightVolume.js

class RCDRenderer extends AbstractRenderer {

    constructor(gl, volume, environmentTexture, options) {
        super(gl, volume, environmentTexture, options);
        Object.assign(this, {
            _light                      : [10, 10, 10],
            _lightDefinitions           : [],
            // _lightType                  : 'distant',
            _stepSize                   : 0.00333,
            _alphaCorrection            : 100,
            _absorptionCoefficient      : 1,
            _scattering                 : 0.5,
            _lightVolumeRatio           : 1,
            _lightToggling              : 0,
            _localSizeX                 : 16,
            _localSizeY                 : 16,
            steps                       : 1,
            _majorant                   : 1,
            _absorptionCoefficientMC    : 16,
            _type                       : 1,
            _rayCastingStepSize         : 0.00333,
            _rayCastingAlphaCorrection  : 100
        }, options);

        this._programs = WebGL.buildPrograms(this._gl, {
            monteCarlo      : SHADERS.RCDMonteCarlo,
            resetPhotons    : SHADERS.RCDResetPhotons,
            rayCasting      : SHADERS.RCDRayCasting,
            diffusion       : SHADERS.RCDDiffusion,
            generate        : SHADERS.RCDGenerate,
            integrate       : SHADERS.RCDIntegrate,
            render          : SHADERS.RCDRender,
            reset           : SHADERS.RCDReset
        }, MIXINS);

        this._transferFunction = WebGL.createTexture(gl, {
            width  : 2,
            height : 1,
            data   : new Uint8Array([255, 0, 0, 0, 255, 0, 0, 255]),
            wrapS  : gl.CLAMP_TO_EDGE,
            wrapT  : gl.CLAMP_TO_EDGE,
            min    : gl.LINEAR,
            mag    : gl.LINEAR
        });

        this._lightDefinitions = [
            new LightDefinition('distant', [10,10,10], true),
            new LightDefinition('point', [220,125,50], false)
        ]

        if (this._volume.ready) {
            this._initVolume();
        }
    }

    _switchToType(type) {
        const gl = this._gl;
        this._type = type;
        if (type === 0) {
            this._initPhotons();
        } else if (this._photonBuffer) {
            gl.deleteBuffer(this._photonBuffer);
            this._photonBuffer = null;
        }
        this._resetLightField();
    }

    _initPhotons() {
        if (!this._lightVolumeDimensions)
            return;
        const gl = this._gl;
        const dimensions = this._lightVolumeDimensions;
        if (this._photonBuffer) {
            gl.deleteBuffer(this._photonBuffer);
        }
        // struct Photon {     //
        //     vec3 position;  // 4 * 4B
        //     vec3 direction; // 4 * 4B
        //     vec3 radiance;  // 4 * 4B
        //     vec3 color;     // 4 * 4B
        //     uint bounces;   // 4B
        //     uint samples;   // 4B
        //          padding    // ??
        // };                  // Where do we get 20 though???
        const bufferSize = 12 * 4 * dimensions.width * dimensions.height * dimensions.depth;
        this._photonBuffer = gl.createBuffer();
        gl.bindBuffer(gl.SHADER_STORAGE_BUFFER, this._photonBuffer);
        gl.bindBufferBase(gl.SHADER_STORAGE_BUFFER, 0, this._photonBuffer);
        gl.bufferData(gl.SHADER_STORAGE_BUFFER, bufferSize, gl.STATIC_DRAW);
    }

    setVolume(volume) {
        this._volume = volume;
        this._initVolume();
        this.reset();
    }

    _initVolume() {
        const volumeDimensions = this._volume.getDimensions('default');
        this._volumeDimensions = volumeDimensions;
        this._setLightVolumeDimensions();
        console.log("Volume Dimensions: " + volumeDimensions.width + " " + volumeDimensions.height + " " + volumeDimensions.depth);
        if (this._type === 0)
            this._initPhotons();
        this._resetLightField();
        this.counter = 0;
    }

    _setLightVolumeDimensions() {
        if (!this._volumeDimensions) {
            return;
        }
        const volumeDimensions = this._volumeDimensions;
        this._lightVolumeDimensions = {
            width: Math.floor(volumeDimensions.width / this._lightVolumeRatio),
            height: Math.floor(volumeDimensions.height / this._lightVolumeRatio),
            depth: Math.floor(volumeDimensions.depth / this._lightVolumeRatio)
        };
        console.log("Light Volume Dimensions: " + this._lightVolumeDimensions.width + " " +
            this._lightVolumeDimensions.height + " " + this._lightVolumeDimensions.depth);
    }

    _resetLightField() {
        if (!this._volumeDimensions) {
            return;
        }
        const gl = this._gl;
        console.log("Reset Light Field")
        if (this._energyDensityVolume) {
            gl.deleteTexture(this._energyDensityVolume);
        }
        this._createLightVolume();
        this._resetDiffusionField();
        switch(this._type) {
            case 0: this.resetPhotons(); break;
            case 1: this._rayCasting(); break;
        }
        this.counter = 0;
    }

    _createLightVolume() {
        const gl = this._gl;
        const dimensions = this._lightVolumeDimensions;
        // Energy density
        this._energyDensityVolume = gl.createTexture();

        // TODO separate function in WebGL.js
        gl.bindTexture(gl.TEXTURE_3D, this._energyDensityVolume);
        gl.texStorage3D(gl.TEXTURE_3D, 1, gl.R32F, dimensions.width, dimensions.height, dimensions.depth);
        gl.texParameteri(gl.TEXTURE_3D, gl.TEXTURE_MIN_FILTER, gl.LINEAR);
        gl.texParameteri(gl.TEXTURE_3D, gl.TEXTURE_MAG_FILTER, gl.LINEAR);
        gl.texParameteri(gl.TEXTURE_3D, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE);
        gl.texParameteri(gl.TEXTURE_3D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);
        gl.texParameteri(gl.TEXTURE_3D, gl.TEXTURE_WRAP_R, gl.CLAMP_TO_EDGE);
    }
    
    _createDiffusionLightVolume() {
        const gl = this._gl;
        const dimensions = this._lightVolumeDimensions;
        // Energy density
        this._energyDensityDiffusion = gl.createTexture();

        // TODO separate function in WebGL.js
        gl.bindTexture(gl.TEXTURE_3D, this._energyDensityDiffusion);
        gl.texStorage3D(gl.TEXTURE_3D, 1, gl.R32F, dimensions.width, dimensions.height, dimensions.depth);
        gl.texParameteri(gl.TEXTURE_3D, gl.TEXTURE_MIN_FILTER, gl.LINEAR);
        gl.texParameteri(gl.TEXTURE_3D, gl.TEXTURE_MAG_FILTER, gl.LINEAR);
        gl.texParameteri(gl.TEXTURE_3D, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE);
        gl.texParameteri(gl.TEXTURE_3D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);
        gl.texParameteri(gl.TEXTURE_3D, gl.TEXTURE_WRAP_R, gl.CLAMP_TO_EDGE);
    }

    _resetDiffusionField() {
        const gl = this._gl;
        console.log("Reset Diffusion Light Field")
        if (this._energyDensityDiffusion)
            gl.deleteTexture(this._energyDensityDiffusion);
        this._createDiffusionLightVolume();
    }

    destroy() {
        const gl = this._gl;
        Object.keys(this._programs).forEach(programName => {
            gl.deleteProgram(this._programs[programName].program);
        });
        if (this._energyDensityVolume)
            gl.deleteTexture(this._energyDensityVolume);
        if (this._energyDensityDiffusion)
            gl.deleteTexture(this._energyDensityDiffusion);
        if (this._photonBuffer)
            gl.deleteBuffer(this._photonBuffer);
        super.destroy();
    }

    _resetFrame() {
        const gl = this._gl;

        const program = this._programs.reset;
        gl.useProgram(program.program);

        gl.drawArrays(gl.TRIANGLE_FAN, 0, 4);
    }

    _monteCarlo() {
        const gl = this._gl;
        const program = this._programs.monteCarlo;
        gl.useProgram(program.program);

        gl.activeTexture(gl.TEXTURE0);
        gl.bindTexture(gl.TEXTURE_3D, this._volume.getTexture());

        gl.activeTexture(gl.TEXTURE1);
        gl.bindTexture(gl.TEXTURE_2D, this._transferFunction);

        gl.uniform1i(program.uniforms.uVolume, 0);
        gl.uniform1i(program.uniforms.uTransferFunction, 1);

        const dimensions = this._lightVolumeDimensions;
        gl.uniform3i(program.uniforms.uSize, dimensions.width, dimensions.height, dimensions.depth);

        gl.uniform1f(program.uniforms.uAbsorptionCoefficient, this._absorptionCoefficientMC)

        gl.uniform1f(program.uniforms.uRatio, Math.floor(this._lightVolumeRatio));

        gl.bindBuffer(gl.SHADER_STORAGE_BUFFER, this._photonBuffer);
        gl.bindBufferBase(gl.SHADER_STORAGE_BUFFER, 0, this._photonBuffer);

        // gl.bindBuffer(gl.SHADER_STORAGE_BUFFER, this._photonBuffer);
        // gl.bindBufferBase(gl.SHADER_STORAGE_BUFFER, 1, this._photonBuffer);

        gl.uniform1f(program.uniforms.uRandSeed, Math.random());
        gl.uniform1f(program.uniforms.uMajorant, this._majorant);
        gl.uniform1ui(program.uniforms.uSteps, this.steps);

        gl.uniform3fv(program.uniforms.light, this._light);

        gl.bindImageTexture(0, this._energyDensityVolume, 0, true, 0, gl.WRITE_ONLY, gl.R32F);

        gl.dispatchCompute(Math.ceil(dimensions.width / this._localSizeX),
            Math.ceil(dimensions.height / this._localSizeY),
            dimensions.depth);
    }

    resetPhotons() {
        const gl = this._gl;
        console.log("Reset Photons")
        const program = this._programs.resetPhotons;
        gl.useProgram(program.program);

        gl.uniform1f(program.uniforms.uRandSeed, Math.random());

        gl.bindBuffer(gl.SHADER_STORAGE_BUFFER, this._photonBuffer);
        gl.bindBufferBase(gl.SHADER_STORAGE_BUFFER, 0, this._photonBuffer);

        const dimensions = this._lightVolumeDimensions;

        gl.uniform3i(program.uniforms.uSize, dimensions.width, dimensions.height, dimensions.depth);
        gl.uniform3fv(program.uniforms.light, this._light);

        gl.dispatchCompute(Math.ceil(dimensions.width / this._localSizeX),
            Math.ceil(dimensions.height / this._localSizeY),
            dimensions.depth);
    }

    _rayCasting() {
        const gl = this._gl;
        const program = this._programs.rayCasting;
        gl.useProgram(program.program);

        gl.activeTexture(gl.TEXTURE0);
        gl.bindTexture(gl.TEXTURE_3D, this._volume.getTexture());

        gl.activeTexture(gl.TEXTURE1);
        gl.bindTexture(gl.TEXTURE_2D, this._transferFunction);

        gl.uniform1i(program.uniforms.uVolume, 0);
        gl.uniform1i(program.uniforms.uTransferFunction, 1);

        const dimensions = this._lightVolumeDimensions;
        gl.uniform3i(program.uniforms.uSize, dimensions.width, dimensions.height, dimensions.depth);

        gl.uniform1f(program.uniforms.uAbsorptionCoefficient, this._absorptionCoefficient);

        gl.uniform3fv(program.uniforms.uLight, this._light);

        gl.uniform1f(program.uniforms.uStepSize, this._rayCastingStepSize);
        gl.uniform1f(program.uniforms.uAlphaCorrection, this._rayCastingAlphaCorrection);

        gl.bindImageTexture(0, this._energyDensityVolume, 0, true, 0, gl.READ_WRITE, gl.R32F);

        gl.dispatchCompute(Math.ceil(dimensions.width / this._localSizeX),
            Math.ceil(dimensions.height / this._localSizeY),
            dimensions.depth);
    }

    _diffusion() {
        const gl = this._gl;

        const program = this._programs.diffusion;
        gl.useProgram(program.program);

        gl.bindImageTexture(0, this._energyDensityVolume, 0, true, 0, gl.READ_ONLY, gl.R32F);
        gl.bindImageTexture(1, this._energyDensityDiffusion, 0, true, 0, gl.READ_WRITE, gl.R32F);

        const dimensions = this._lightVolumeDimensions;

        gl.uniform3i(program.uniforms.uSize, dimensions.width, dimensions.height, dimensions.depth);
        gl.uniform1f(program.uniforms.scattering, this._scattering);

        gl.uniform1f(program.uniforms.uRatio, Math.floor(this._lightVolumeRatio));

        gl.dispatchCompute(Math.ceil(dimensions.width / this._localSizeX),
            Math.ceil(dimensions.height / this._localSizeY),
            dimensions.depth);
    }

    _generateFrame() {
        const gl = this._gl;
        if (this._type === 0) {
            this._monteCarlo();
        }
        this._diffusion();

        const program = this._programs.generate;
        gl.useProgram(program.program);

        gl.activeTexture(gl.TEXTURE0);
        gl.bindTexture(gl.TEXTURE_3D, this._volume.getTexture());
        gl.activeTexture(gl.TEXTURE1);
        gl.bindTexture(gl.TEXTURE_2D, this._transferFunction);
        gl.activeTexture(gl.TEXTURE2);
        gl.bindTexture(gl.TEXTURE_3D, this._energyDensityVolume);
        gl.activeTexture(gl.TEXTURE3);
        gl.bindTexture(gl.TEXTURE_3D, this._energyDensityDiffusion);

        gl.uniform1i(program.uniforms.uVolume, 0);
        gl.uniform1i(program.uniforms.uTransferFunction, 1);
        gl.uniform1i(program.uniforms.uEnergyDensity, 2);
        gl.uniform1i(program.uniforms.uEnergyDensityDiffusion, 3);
        gl.uniform1f(program.uniforms.uStepSize, this._stepSize);
        gl.uniform1f(program.uniforms.uAlphaCorrection, this._alphaCorrection);
        gl.uniform1f(program.uniforms.uOffset, Math.random());
        gl.uniformMatrix4fv(program.uniforms.uMvpInverseMatrix, false, this._mvpInverseMatrix.m);

        gl.drawArrays(gl.TRIANGLE_FAN, 0, 4);
    }

    _integrateFrame() {
        const gl = this._gl;

        const program = this._programs.integrate;
        gl.useProgram(program.program);

        gl.activeTexture(gl.TEXTURE0);
        gl.bindTexture(gl.TEXTURE_2D, this._accumulationBuffer.getAttachments().color[0]);
        gl.activeTexture(gl.TEXTURE1);
        gl.bindTexture(gl.TEXTURE_2D, this._frameBuffer.getAttachments().color[0]);

        gl.uniform1i(program.uniforms.uAccumulator, 0);
        gl.uniform1i(program.uniforms.uFrame, 1);

        gl.drawArrays(gl.TRIANGLE_FAN, 0, 4);
    }

    _renderFrame() {
        const gl = this._gl;

        const program = this._programs.render;
        gl.useProgram(program.program);

        gl.activeTexture(gl.TEXTURE0);
        gl.bindTexture(gl.TEXTURE_2D, this._accumulationBuffer.getAttachments().color[0]);

        gl.uniform1i(program.uniforms.uAccumulator, 0);

        gl.drawArrays(gl.TRIANGLE_FAN, 0, 4);
    }

    _getProgramFromLightType(index) {
        switch(this._lightVolumes[index].getType()) {
            case 'distant': return this._programs.convection;
            case 'point': return this._programs.convectionPL;
        }
    }

    _getFrameBufferSpec() {
        const gl = this._gl;
        return [{
            width          : this._bufferSize,
            height         : this._bufferSize,
            min            : gl.NEAREST,
            mag            : gl.NEAREST,
            format         : gl.RGBA,
            internalFormat : gl.RGBA,
            type           : gl.UNSIGNED_BYTE
        }];
    }

    _getAccumulationBufferSpec() {
        const gl = this._gl;
        return [{
            width          : this._bufferSize,
            height         : this._bufferSize,
            min            : gl.NEAREST,
            mag            : gl.NEAREST,
            format         : gl.RGBA,
            internalFormat : gl.RGBA,
            type           : gl.UNSIGNED_BYTE
        }];
    }

}