// #package js/main
// #include ../Pixel.js

class KMeansImplementation {

    static weightForMagnitude = 0.8;
    static maxIterations = 1000;

    static findPeaks(data, numOfPeaks) {
        let tempData = data;
        let k = numOfPeaks;
        let centers = this.determineStartingPoints(tempData, k);
        
        let nextIter = true;
        let iterations = 0;
        // repeat until no change anymore
        while (nextIter) {
            iterations++;
            tempData = this.assignToClosestCenter(tempData, centers);
            let newCenters = this.determineNewCenters(tempData, centers);
            centers = newCenters.centers;
            console.log(newCenters);
            if (!newCenters.change || iterations === this.maxIterations)
                nextIter = false;
        }
        console.log('it took ' + iterations + ' iterations');
        console.log(centers);
        return centers;
    }

    // determine starting k points
    static determineStartingPoints(data, k) {
        let points = [];
        points.push(this.getDataAtIndex(data, 0, 1));
        while (points.length != k) {
            let x = Math.floor(Math.random() * (256));
            let y = Math.floor(Math.random() * (256));
            //let furthestPixel = {
            //    pixel: null,
            //    distance: -1
            //};
            //data.forEach(pixel => {
            //    if (pixel.grayscale == 0)
            //        return;
            //    let distance = 0;
            //    points.forEach(p => {
            //        distance += this.calculateDistance(pixel, p);
            //    });
            //    if (furthestPixel.distance < distance || furthestPixel.distance == -1) {
            //        furthestPixel.pixel = pixel;
            //        furthestPixel.distance = distance;
            //    }
            //});
            let pixel = this.getDataAtIndex(data, x, y);
            if (pixel.grayscale != 0) {
                points.push(pixel);
            }
            //points.push(furthestPixel.pixel);
        }
        console.log('starting points are:');
        console.log(points);
        return points;
    }

    static assignToClosestCenter(data, centers) {
        data.forEach(pixel => {
            if (pixel.grayscale === 0)
                return;
            let closestCenter = {
                weight: -1, 
                center: null
            };
            centers.forEach(c => {
                let weight = this.calculateWeight(c, pixel);
                if (weight < closestCenter.weight || closestCenter.weight == -1) {
                    closestCenter = {
                        center: c,
                        weight: weight
                    };
                }
            });
            pixel.center = closestCenter.center;
        });
        return data;
    }

    // calculate centroid of new centers
    static determineNewCenters(data, centers) {
        let change = false;

        let values = [];
        let newCenters = [];
        centers.forEach(c => {
            values.push({
                center: c,
                number: 0,
                x: 0,
                y: 0
            });
        });

        data.forEach(p => {
            if (p.grayscale === 0)
                return;
            values.forEach(v => {
                if (v.center.equals(p.center)) {
                    v.number++;
                    v.x += p.x;
                    v.y += p.y;
                }
            });
        });

        values.forEach(v => {
            let newX = Math.round(v.x / v.number);
            let newY = Math.round(v.y / v.number);
            newCenters.push(this.getDataAtIndex(data, newX, newY));
        });

        // detect any change
        newCenters.forEach(nc => {
            if (centers.every(c => !nc.equals(c)))
                change = true;
        })

        return {
            centers: newCenters,
            change: change
        };
    }

    static calculateWeight(p1, p2) {
        let magnitudePart = this.weightForMagnitude * Math.pow(p1.grayscale - p2.grayscale, 2);
        let distancePart = (1 - this.weightForMagnitude) * this.calculateDistance(p1, p2);
        return magnitudePart + distancePart;
    }

    static calculateDistance(p1, p2) {
        return Math.sqrt(Math.pow(p1.x - p2.x, 2) + Math.pow(p1.y - p2.y, 2));
    }

    static getDataAtIndex(data, x, y) {
        return data[y * 256 + x];
    }

}