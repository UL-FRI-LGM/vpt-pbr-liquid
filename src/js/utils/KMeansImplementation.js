// #package js/main
// #include ../Pixel.js

class KMeansImplementation {

    static weightForMagnitude = 0.6;
    static maxIterations = 1000;
    // total number of pixels = 256 * 256 = 65.536
    static minClusterDistance = 5;
    static minClusters = 3;
    static maxClusters = 8;
    static minClusterSize = 20;

    static findPeaks(data, numOfPeaks) {
        let tempData = data;
        let k = numOfPeaks;
        let centers = this.determineStartingClusterCenters(tempData, k);
        
        let nextIter = true;
        let iterations = 0;
        // repeat until no change anymore
        while (nextIter) {
            iterations++;
            tempData = this.assignToClosestCenter(tempData, centers);
            let newCenters = this.determineNewCenters(tempData, centers);
            centers = newCenters.centers;
            if (!newCenters.change || iterations === this.maxIterations)
                nextIter = false;
        }
        console.log('it took ' + iterations + ' iterations');
        console.log(centers);
        return centers;
    }

    // determine starting k points
    static determineStartingClusterCenters(data, k) {
        let points = [];
        // select first pixel that value of grayscale is not equal 0 and has no value of 0, 254 or 255
        for (const pixel of data) {
                if (pixel.x === 0 || pixel.x === 254 || pixel.x === 255 || pixel.grayscale === 0) {
                    continue;
                }
                points.push(pixel);
                break;
        }
        while (points.length != k) {
            let furthestPixel = {
                pixel: null,
                distance: -1
            };
            data.forEach(pixel => {
                if (pixel.grayscale === 0 || pixel.x === 0 || pixel.x === 254 || pixel.x === 255 || points.includes(pixel))
                    return;
                let distance = 0;
                points.forEach(p => {
                    distance += this.calculateDistance(pixel, p);
                });
                if (furthestPixel.distance < distance || furthestPixel.distance === -1) {
                    furthestPixel.pixel = pixel;
                    furthestPixel.distance = distance;
                }
            });
            points.push(furthestPixel.pixel);
        }
        return points;
    }

    static assignToClosestCenter(data, centers) {
        data.forEach(pixel => {
            if (pixel.grayscale === 0 || pixel.x === 0 || pixel.x === 254 || pixel.x === 255)
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

        let centerInfo = [];
        let newCenters = [];
        // for every center initialize an object that stores information about how many points are assigned to that center 
        // and the accumulate x and y coordinates of them to calculate centroids
        centers.forEach(c => {
            centerInfo.push({
                center: c,
                number: 1,
                x: c.x,
                y: c.y
            });
        });

        data.forEach(p => {
            if (p.grayscale === 0 || p.x === 0 || p.x === 254 || p.x === 255) {
                return;
            }
            // for every pixel check in which center does it belong and add that pixel's data
            centerInfo.forEach(c => {
                if (c.center.equals(p.center)) {
                    c.number++;
                    c.x += p.x;
                    c.y += p.y;
                }
            });
        });

        // calculate new centers as average of all associated pixels
        centerInfo.forEach(c => {
            let newX = Math.round(c.x / c.number);
            let newY = Math.round(c.y / c.number);
            newCenters.push({
                center: this.getDataAtIndex(data, newX, newY),
                number: c.number
            });
        });
        let adjusted = this.adjustClusters(data, newCenters);

        // detect any change
        adjusted.forEach(nc => {
            if (centers.every(c => !nc.center.equals(c)))
                change = true;
        })

        return {
            centers: adjusted.map(nc => nc.center),
            change: change
        };
    }

    // ISODATA extention
    // centers = list of new centers, containing cluster center pixel and number of elements it has
    static adjustClusters(data, centers) {
        // check if too many/few clusters
        if (centers.length > this.maxClusters) {
            let nearestClusterCenters = this.determineNearestClusters(centers);
            centers = this.mergeClusters(data, centers, nearestClusterCenters.center1, nearestClusterCenters.center2);
        }
        if (centers.length < this.minClusters) {
            let furthestClusterCenter = this.determineFurthestCluster(centers);
            centers = this.splitCluster(data, centers, furthestClusterCenter);
        }
        // check if cluster centers are too close to each other, than merge
        for (const c1 of centers) {
            let nearestCluster = this.determineNearestCluster(centers, c1);
            if (c1.number < this.minClusterSize) {
                centers = this.mergeClusters(data, centers, nearestCluster.center1, nearestCluster.center2);
            }
            if (nearestCluster.distance < this.minClusterDistance) {
                centers = this.mergeClusters(data, centers, nearestCluster.center1, nearestCluster.center2);
            }
        }        
        return centers;
    }

    static determineNearestCluster(centers, clusterCenter) {
        let nearestCluster = {
            center1: clusterCenter,
            center2: null,
            distance: -1
        };
        centers.forEach(c => {
            if (c.center.equals(clusterCenter)) {
                return;
            }
            let distance = this.calculateDistance(c.center, clusterCenter);
            if (nearestCluster.distance === -1 || distance < nearestCluster.distance) {
                nearestCluster.center2 = c;
                nearestCluster.distance = distance;
            }
        });
        return nearestCluster;
    }

    static determineNearestClusters(centers) {
        let nearestClusters = {
            center1: null,
            center2: null,
            distance: -1
        }
        centers.forEach(c1 => {
            let nearestCluster = this.determineNearestCluster(centers, c1);
            if (nearestClusters.distance === -1 || nearestCluster.distance < nearestClusters.distance) {
                nearestClusters.center1 = nearestCluster.center1;
                nearestClusters.center2 = nearestCluster.center2;
                nearestClusters.distance = nearestCluster.distance;
            }
        });
        return nearestClusters;
    }

    static determineFurthestCluster(centers) {
        let longest = {
            center: null,
            distance: -1
        }
        centers.forEach(c1 => {
            let totalDistance = 0;
            centers.forEach(c2 => {
                if (c1.center.equals(c2.center))
                    return;
                totalDistance += this.calculateDistance(c1, c2)
            });
            if (totalDistance > longest.distance || longest.distance === -1) {
                longest.center = c1;
                longest.distance = totalDistance;
            }
        });
        return longest.center;
    }

    // remove c2, assign everything to c1
    static mergeClusters(data, centers, c1, c2) {
        data.forEach(p => {
            if (p.center !== null && p.center.equals(c2.center)) {
                p.center = c1.center;
            }
        });
        return centers.filter(c => !c.center.equals(c2.center));
    }

    // split cluster by y axis on its center
    static splitCluster(data, centers, ce) {
        let c1 = {
            center: null,
            number: 0,
            x: 0,
            y: 0
        };
        let c2 = {
            center: null,
            number: 0,
            x: 0,
            y: 0
        };
        data.forEach(p => {
            if (p.center !== null && p.center.equals(ce.center)) {
                // check in which new cluster it belongs
                if (p.y >= ce.center.y) {
                    p.center = c1;
                    c1.number++;
                    c1.x += p.x;
                    c1.y += p.y;
                } else {
                    p.center = c2;
                    c2.number++;
                    c2.x += p.x;
                    c2.y += p.y;
                }
            }
        });
        c1.center = this.determineNewCenter(data, c1).center;
        c2.center = this.determineNewCenter(data, c2).center;
        let newCenters = centers.filter(c => !c.center.equals(ce.center));
        newCenters.push(c1, c2);
        return newCenters;
    }

    static determineNewCenter(data, c) {
        let newX = Math.round(c.x / c.number);
        let newY = Math.round(c.y / c.number);
        return {
            center: this.getDataAtIndex(data, newX, newY),
            number: c.number
        }
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