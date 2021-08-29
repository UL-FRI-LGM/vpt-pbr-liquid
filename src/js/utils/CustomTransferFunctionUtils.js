// #package js/main
// #include ../Pixel.js
// #include Colormap.js
// #include KMeansImplementation.js

class CustomTransferFunctionUtils {

    static createCustomTf(tfArray) {
        let grayscaledPixels = this.convertToGrayscale(tfArray);

        let coloredPixels = this.findAndDetectPeaks(grayscaledPixels);

        let colors = coloredPixels.map(pixel => pixel.color);

        let customTf = new Array(colors.length * 4);
        for (let i = 0; i < colors.length; i++) {
            customTf[4 * i] = colors[i][0];
            customTf[4 * i + 1] = colors[i][1];
            customTf[4 * i + 2] = colors[i][2];
            customTf[4 * i + 3] = colors[i][3];
        }

        return customTf;
    }

    static convertToGrayscale(pixelData) {
        let pixels = [];
        pixelData.forEach((element, i, array) => {
            if (i % 4 !== 0)
                return;
            let index = i / 4;
            let x = index % 256;
            let y = parseInt(index / 256);
            let r = this.changeInfinityToMax(array[i], 255);
            let g = this.changeInfinityToMax(array[i + 1], 255);
            let b = this.changeInfinityToMax(array[i + 2], 255);
            let a = this.changeInfinityToMax(array[i + 3], 255);

            let grayscale = 255 - (r + g + b) / 3;
            pixels.push(new Pixel(x, y, r, g, b, a, grayscale));
        });
        return pixels;
    }

    static changeInfinityToMax(value, newValue) {
    if (value === Number.POSITIVE_INFINITY)
            return newValue;
        return value;
    }

    static findAndDetectPeaks(data) {
        let peaks = this.findPeaks(data);
        console.log('found peaks:');
        console.log(peaks.length);
        this.colorPeaksAndPixels(data, peaks, 0.1, 0.9);
        return data;
    }

    static findPeaks(data) {
        return KMeansImplementation.findPeaks(data, 5);
    }

    static getDataAtIndex(data, x, y) {
        return data[y * 256 + x];
    }

    static colorPeaksAndPixels(data, peaks, start, end) {
        let colormap = Colormap.returnColormap('blue');
        let interval = (end - start) / peaks.length;
        peaks.forEach((pixel, i) => {
            pixel.color = this.generatePeakColor(colormap, interval, start, i);
        });
        data.forEach(pixel => {
            if (pixel.peak)
                return;
            else if (pixel.x === 254)
                // cube (green)
                pixel.color = [0, 255, 0, 255];
            else if (pixel.x === 255) 
                // floor (yellow)
                pixel.color = [255, 255, 0, 255];
            else if (pixel.x === 0 || pixel.grayscale === 0)
                // air (black with no alpha)
                pixel.color = [255, 255, 255, 0];
            else {
                // not floor, cube or peak, interpolate
                let color = ColorInterpolationUtils.getColorRatios(pixel, peaks);
                pixel.color = color;
            }
        });
    }

    static generatePeakColor(colormap, interval, start, i) {
        let length = colormap.length;
        let positionInInterval = Math.random() * interval;
        let positionOnScale = i * interval + positionInInterval;
        let positionInArray = Math.round(positionOnScale * length);
        let color = colormap[positionInArray];
        return [Math.round(color[0] * 255), Math.round(color[1] * 255), Math.round(color[2] * 255)];
    }

}