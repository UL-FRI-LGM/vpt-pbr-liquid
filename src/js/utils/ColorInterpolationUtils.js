// #package js/main

class ColorInterpolationUtils {

    static getColorRatios(pixel, peaks) {
        colorRatios = Array(peaks.length).fill(1);
        peaks.forEach((peak1, i1) => {
            peaks.forEach((peak2, i2) => {
                if (i1 == i2)
                    return;
                d = this.calculateProjectionDistance(peak1, peak2, pixel);
                colorRatios[index1] *= this.limit(d);
            });
        });
        totalRatiosSum = 0;
        colorRatios.forEach(ratio => totalRatiosSum += ratio);
        colorRatios.forEach(ratio => ratio /= totalRatiosSum);
        color = this.getColorMix(peaks, colorRatios);
        return color;
    }

    static calculateProjectionDistance(p1, p2, pixel) {
        let k2 = p2.x * p2.x - p2.x * p1.x + p2.y * p2.y - p2.y * p1.y;
        let k1 = p1.x * p1.x - p2.x * p1.x + p1.y * p1.y - p2.y * p1.y;
        let ab2 = (p1.x - p2.x) * (p1.x - p2.x) + (p1.y - p2.y) * (p1.y - p2.y);
        let kcom = pixel.x * (p1.x - p2.x) + pixel.y * (p1.y - p2.y);
        let d1 = (k1 - kcom) / ab2;
        let d2 = (k2 + kcom) / ab2;
        return d2;
    }

    static limit(value) {
        if (value < 0)
            return 0;
        if (value > 1)
            return 1;
        return value;
    }

    static getColorMix(peaks, colorRatios) {
        let r, g, b = 0;
        peaks.forEach((peak, index) => {
            r += peak.color[0] * colorRatios[index];
            g += peak.color[1] * colorRatios[index];
            b += peak.color[2] * colorRatios[index];
        });
        return [Math.round(r), Math.round(g), Math.round(b), 255];
    }

} 