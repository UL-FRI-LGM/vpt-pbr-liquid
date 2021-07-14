// #package js/main

class Pixel {

    constructor(x, y, r, g, b, a, grayscale) {
        this.x = x;
        this.y = y;
        this.color = [r, g, b, a];
        this.peak = false;
        this.grayscale = grayscale;
    }

    getColor() {
        return this.color;
    }

}