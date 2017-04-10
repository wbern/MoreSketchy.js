/* eslint-disable no-param-reassign, no-unused-vars */
const munkres = require('munkres-js');

// ////////////////////////////////////////
//              Sketchy.js
// ////////////////////////////////////////
//
// JavaScript Shape Matching Works best with Raphael SketchPad Development
// started March 2013
//

// Immediately invoke an anonymous function to keep the global scope clean.
// Parameters:
// - global: will be the global object; called with "this" from global scope
// - undefined: keeps "undefined" undefined; no 2nd arg will make it undefined
// Namespace everything
const Sketchy = {};

/* Jordan's Algorithms */
// Test function for front-end application development
Sketchy.randomShapeMatch = function (shape1, shape2) {
  return Math.random();
};

/* Kyle's Algorithms */
// Takes in SVG data (from Raphael SketchPad) and outputs an array of paths,
// each of which is an array of points in {x: Number, y: Number} format. This is
// useful for preprocessing for Simplify.js or any other algorithm operating on
// simple paths.
Sketchy.convertSVGtoPointArrays = function (json) {
  let splitPath;
  let point;
  const paths = [];

  json = JSON.parse(json);
  for (let i = 0; i < json.length; i++) {
    // Take the SVG data for the current path, cut off the M at the beginning, and
    // then explode the string into an array, split at the "L" character.  This is
    // the format from Raphael SketchPad
    splitPath = null;
    if (typeof json[i] === 'string') {
      splitPath = json[i]
        .slice(1)
        .split(/[A-z]/);
    } else {
      splitPath = json[i]
        .path
        .slice(1)
        .split(/[A-z]/);
    }
    paths[i] = [];
    for (let j = 0; j < splitPath.length; j++) {
      point = splitPath[j].split(',');
      paths[i][j] = {
        x: parseInt(point[0], 10),
        y: parseInt(point[1], 10)
      };
    }
  }
  return paths;
};
// Takes in an array of paths, each of which is an array of points in {x:
// Number, y: Number} format, and outputs it in Raphael SketchPad-style JSON/SVG
// data.  Essentially reverses the above and makes the same drawing decisions as
// Raphael SketchPad (e.g. black, stroke-width of 5).
Sketchy.convertPointArraysToSVG = function (paths) {
  const json = [];

  for (let i = 0; i < paths.length; i++) {
    json[i] = {
      fill: 'none',
      stroke: '#000000',
      path: 'M',
      'stroke-opacity': 1,
      'stroke-width': 5,
      'stroke-linecap': 'round',
      'stroke-linejoin': 'round',
      transform: [],
      type: 'path'
    };
    json[i].path += `${paths[i][0].x},${paths[i][0].y}`;
    for (let j = 1; j < paths[i].length; j++) {
      json[i].path += `L${paths[i][j].x},${paths[i][j].y}`;
    }
  }
  return JSON.stringify(json); // TODO: better distinguish between JSON strings and objects
};

// Takes in SVG data (from Raphael SketchPad) and outputs an svgXML file.
Sketchy.convertSVGtoXML = function (jsonParam, parsed) {
  let splitPath;
  let point;
  let svgXML;
  let json;
  // svgXML = "<?xml version=\"1.0\" standalone=\"no\"?>\n<!DOCTYPE svg PUBLIC
  // \"-//W3C//DTD SVG 1.1//EN\"
  // \"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">\n";
  svgXML = '<svg>'; // width=\"100%\" height=\"100%\" version=\"1.1\" xmlns=\"http://www.w3.org/2000/svg\"> ";
  if (!parsed) {
    json = JSON.parse(jsonParam);
  } else {
    json = jsonParam;
  }
  for (let i = 0; i < jsonParam.length; i++) {
    svgXML += '\n<path fill="none" stroke-opacity="1" stroke="#000000" stroke-linecap="round" s' +
        'troke-width="5" stroke-linejoin="round" type="path" d="M ';
    splitPath = jsonParam[i]
      .path
      .slice(1)
      .split(/(?=[L])/);
    for (let j = 0; j < splitPath.length; j++) {
      point = splitPath[j].split(',');
      svgXML += `${point[0]} ${point[1]} `;
    }
    svgXML += '"/>';
  }
  svgXML += '\n</svg>';
  return svgXML;
};

// shape1 and shape2 should be stringified JSON data from Raphael SketchPad
Sketchy.shapeContextMatch = function (shape1, shape2, convertToPoints = true, accuracyPercentage = 10) {
  const pointsPerShape = parseInt(accuracyPercentage / 2, 10); // constant... 25 is pretty fast... 50 is probably best
  let points1;
  let points2;

    // 0.125 gives a bin out to points 2x the average
  let distanceBinSmallest = 0.125;
  const distanceBinCount = 5;

  const angleBinCount = 12;

  // Scatter points around each of the paths.  The algorithm will only be using
  // these points (as feature descriptors), not the shapes.
  if (convertToPoints) {
    points1 = Sketchy.scatterPoints(Sketchy.convertSVGtoPointArrays(shape1), pointsPerShape);
    points2 = Sketchy.scatterPoints(Sketchy.convertSVGtoPointArrays(shape2), pointsPerShape);
  } else {
    points1 = Sketchy.scatterPoints(shape1, pointsPerShape);
    points2 = Sketchy.scatterPoints(shape2, pointsPerShape);
  }

  const distanceMatrix1 = [];
  const distanceMatrix2 = [];

  // Create a square 2D array and initialize it with 0s in the diagonal
  for (let i = 0; i < pointsPerShape; i++) {
    distanceMatrix1[i] = [];
    distanceMatrix1[i][i] = 0;
    distanceMatrix2[i] = [];
    distanceMatrix2[i][i] = 0;
  }

  // Go through the upper triangle of the matrix, computing the distance,
  // mirroring to the lower
  let distanceTotal1 = 0;
  let distanceTotal2 = 0;
  for (let i = 0; i < pointsPerShape - 1; i++) {
    for (let j = i + 1; j < pointsPerShape; j++) {
      distanceMatrix1[i][j] = Sketchy.euclideanDistance(points1[i].x, points1[i].y, points1[j].x, points1[j].y);
      distanceMatrix1[j][i] = distanceMatrix1[i][j];
      distanceMatrix2[i][j] = Sketchy.euclideanDistance(points2[i].x, points2[i].y, points2[j].x, points2[j].y);
      distanceMatrix2[j][i] = distanceMatrix2[i][j];
      distanceTotal1 += distanceMatrix1[i][j];
      distanceTotal2 += distanceMatrix2[i][j];
    }
  }
  distanceTotal1 *= 2; // 0s already summed in, we just need to double it since we only went through the upper triangle
  distanceTotal2 *= 2;
  const distanceMean1 = distanceTotal1 / (pointsPerShape ** 2);
  const distanceMean2 = distanceTotal2 / (pointsPerShape ** 2);

  // Normalize by the mean distance.  This achieves scale invariance. Translation
  // invariance is inherent to the fact that distance measurements are made
  // relative to each point.
  for (let i = 0; i < pointsPerShape; i++) {
    for (let j = 0; j < pointsPerShape; j++) {
      distanceMatrix1[i][j] /= distanceMean1;
      distanceMatrix2[i][j] /= distanceMean2;
    }
  }

  const distanceBins1 = [];
  const distanceBins2 = [];

  // Initialize the distance bins with all 0s
  for (let i = 0; i < pointsPerShape; i++) {
    distanceBins1[i] = [];
    distanceBins2[i] = [];
    for (let j = 0; j < pointsPerShape; j++) {
      distanceBins1[i][j] = 0;
      distanceBins2[i][j] = 0;
    }
  }

  // Double the acceptable radius each iteration, increasing the bin number each
  // time a point is still in the running.  0 means the point was not in any bins
  // (and will not be counted), 1 means it was in the outer, and distanceBinCount
  // (e.g. 5) means it is in the closest bin (including the same point).
  for (let k = 0; k < distanceBinCount; k++) {
    for (let i = 0; i < pointsPerShape; i++) {
      for (let j = i + 1; j < pointsPerShape; j++) {
        if (distanceMatrix1[i][j] < distanceBinSmallest) {
          distanceBins1[i][j]++;
          distanceBins1[j][i]++;
        }
        if (distanceMatrix2[i][j] < distanceBinSmallest) {
          distanceBins2[i][j]++;
          distanceBins2[j][i]++;
        }
      }
    }
    distanceBinSmallest *= 2;
  }

  const angleMatrix1 = [];
  const angleMatrix2 = [];

  // Angles // Create a square 2D array and initialize it with 0s in the diagonal
  for (let i = 0; i < pointsPerShape; i++) {
    angleMatrix1[i] = [];
    angleMatrix2[i] = [];
    angleMatrix1[i][i] = 0;
    angleMatrix2[i][i] = 0;
  }

  // Compute the angle matrix, much like the distance matrix
  for (let i = 0; i < pointsPerShape - 1; i++) {
    for (let j = i + 1; j < pointsPerShape; j++) {
      // Adding 2pi and modding by 2pi changes the -pi to pi range to a 0 to 2pi range
      angleMatrix1[i][j] = (
        Math.atan2(points1[j].y - points1[i].y, points1[j].x - points1[i].x) + (2 * Math.PI)) % ((2 * Math.PI)
      );
      angleMatrix2[i][j] = (
        Math.atan2(points2[j].y - points2[i].y, points2[j].x - points2[i].x) + (2 * Math.PI)) % ((2 * Math.PI)
      );

      // The matrix is somewhat mirrored over the diagonal, but angles must be flipped
      // around
      angleMatrix1[j][i] = (angleMatrix1[i][j] + Math.PI) % (2 * Math.PI);
      angleMatrix2[j][i] = (angleMatrix2[i][j] + Math.PI) % (2 * Math.PI);
    }
  }

  // Initialize the angle bins
  const angleBins1 = [];
  const angleBins2 = [];

  for (let i = 0; i < pointsPerShape; i++) {
    angleBins1[i] = [];
    angleBins2[i] = [];
  }

  // Compute the angle bins
  // TODO: save efficiency by automatically calculating mirror by adding
  // angleBinCount/2 then modding by angleBinCount?
  for (let i = 0; i < pointsPerShape; i++) {
    for (let j = 0; j < pointsPerShape; j++) {
      angleBins1[i][j] = 1 + Math.floor(angleMatrix1[i][j] / ((2 * Math.PI) / angleBinCount));
      angleBins2[i][j] = 1 + Math.floor(angleMatrix2[i][j] / ((2 * Math.PI) / angleBinCount));
    }
  }

  // Cost Matrix // Compute the cost matrix.  This skips the combined histogram
  // for the sake of efficiency.  TODO: make more efficient by only calculating
  // upper triangle
  const costMatrix = [];

  for (let i = 0; i < pointsPerShape; i++) {
    costMatrix[i] = [];
    let ksum;
    let compare;
    let logr;
    let theta;

    for (let j = 0; j < pointsPerShape; j++) {
      // Go through all K bins.
      ksum = 0;

      for (logr = 1; logr <= distanceBinCount; logr++) {
        for (theta = 1; theta <= angleBinCount; theta++) {
          // calculate hik and hjk
          const hik = Sketchy.shapeContextHistogram(i, logr, theta, distanceBins1, angleBins1);
          const hjk = Sketchy.shapeContextHistogram(j, logr, theta, distanceBins2, angleBins2);
          compare = (hik + hjk === 0)
            ? 0
            : (((hik - hjk) ** 2) / (hik + hjk));
          ksum += compare;
        }
      }
      costMatrix[i][j] = (1 / 2) * ksum;
    }
  }

  // Normalize total cost by the number of points per shape.
  // let result = Sketchy.hungarian(costMatrix, false, true) / pointsPerShape;
  let result = munkres(costMatrix) / pointsPerShape; // Rely on maintained dependency

  // Convert total error to a percentage Modify the constant below (originally
  // 0.175) to modify how sensitive this function is to error.  Higher numbers
  // make it more forgiving.
  // Note: this is Gaussian function.
  result = Math.exp(-((result * result) / 0.175));

  return result;
};

// Sums up the number of points (relative to point pointIndex) in a particular
// bin, defined by distanceBinNumber and angleBinNumber.  The pair
// (distanceBinNumber, angleBinNumber) defines what is typically called k, the
// polar bin.  This replaces the space requirement of a 2D/k-bin histogram for
// each point.
Sketchy.shapeContextHistogram = function (pointIndex, distanceBinNumber, angleBinNumber, distanceBins, angleBins) {
  let accumulator = 0;
  const numberOfPoints = distanceBins.length;
  for (let i = 0; i < numberOfPoints; i++) {
    if (i !== pointIndex &&
    distanceBins[pointIndex][i] === distanceBinNumber &&
    angleBins[pointIndex][i] === angleBinNumber) {
      accumulator++;
    }
  }
  // Normalize by numberOfPoints (technically should be by numberOfPoints-1?)
  // Shouldn't make a difference
  return accumulator / numberOfPoints;
};

// Compute the Euclidean distance (as a crow flies) between two points. Shortest
// distance between two pixels
Sketchy.euclideanDistance = function (x1, y1, x2, y2) {
  return Math.sqrt(((x1 - x2) ** 2) + ((y1 - y2) ** 2));
};

// Compute the city block distance (or Manhattan distance) between two points.
// Shortest 4-connected path between two pixels
Sketchy.cityBlockDistance = function (x1, y1, x2, y2) {
  return Math.abs(x1 - x2) + Math.abs(y1 - y2);
};

// Compute the chessboard distance between two points.
Sketchy.chessboardDistance = function (x1, y1, x2, y2) {
  return Math.max(Math.abs(x1 - x2), Math.abs(y1 - y2));
};

// Compute the length of a path. Given an array of points in {x: Number, y:
// Number} format, calculate the sum of the distances between consecutive
// points.  The distance function must be specified.
// TODO: Currently, there is no error checking (e.g. a valid callback). Either
// add it or make private.
Sketchy.computeLength = function (path, distanceFunction) {
  let distance = 0;
  for (let i = 0; i < path.length - 1; i++) {
    distance += distanceFunction(path[i].x, path[i].y, path[i + 1].x, path[i + 1].y);
  }
  return distance;
};

// Accepts an array of point arrays (multiple paths) and distributes a specified
// number of points accross them, using the distributePointsAcrossPath method.
// This returns numberOfPoints point objects in a single array, thus, path
// information is intentionally lost.
Sketchy.scatterPoints = function (paths, numberOfPoints) {
  let pointsNotAssigned = numberOfPoints;
  const result = [];
  let pathLength;
  let lengthNotCovered;
  let numberOfPointsForPath;
  let path;
  let point;

  // Compute the length of all paths
  lengthNotCovered = 0;
  for (let i = 0; i < paths.length; i++) {
    lengthNotCovered += Sketchy.computeLength(paths[i], Sketchy.euclideanDistance);
  }

  // Scatter points
  for (let i = 0; i < paths.length; i++) {
    path = paths[i];

    // Determine how many points this path will get, based on distance The last path
    // automatically gets any remaining points just in case there is imprecision in
    // the calculations
    pathLength = Sketchy.computeLength(path, Sketchy.euclideanDistance);
    numberOfPointsForPath = Math.round((pathLength / lengthNotCovered) * pointsNotAssigned);
    if (i === paths.length - 1) {
      path = Sketchy.distributePointsAcrossPath(path, pointsNotAssigned);
      pointsNotAssigned = 0;
      lengthNotCovered = 0;
    } else {
      path = Sketchy.distributePointsAcrossPath(path, numberOfPointsForPath);
      pointsNotAssigned -= numberOfPointsForPath;
      lengthNotCovered -= pathLength;
    }

    // Put the points into the result array, disregarding separate paths
    for (let j = 0; j < path.length; j++) {
      point = path[j];
      result.push({ x: point.x, y: point.y }); // copy of the point, not reference
    }
  }

  return result;
};

// Old version of algorithm that selects points from the original list of points
// at a fixed interval... [1,1,1,2,3,4,7,8,9]  <--- select 5 points      ^   ^
// ^   ^   ^ This is not ideal because points tend to get bunched up.  Use the
// below, uncommented implementation instead. Sketchy.distributePointsAcrossPath
// = function(path, numberOfPoints) {   var result, pathIndexDelta, point, i,
// currPathIndex=0;   if(numberOfPoints <= 0) {     return [];   }
// if(numberOfPoints === 1) {     point = path[Math.floor((path.length-1)/2)];
// // reference to original     return [{x:point.x, y:point.y}]; // return a
// copy   }   pathIndexDelta = path.length/(numberOfPoints-1);   // If
// numberOfPoints >= 2, we will manually add the first and last points   // Add
// the first   point = path[0];   result = [{x:point.x, y:point.y}];   for(i=1;
// i<numberOfPoints-1; i++) { currPathIndex += pathIndexDelta;     point =
// path[Math.round(currPathIndex)];     result.push({x:point.x, y:point.y}); //
// TODO: an error occurs (point is undefined) here when a short paths are drawn
// and shapeContextMatch is called }   // Add the last   point =
// path[path.length-1];   result.push({x:point.x, y:point.y});   return result;
// }; Turn an array of points into another array representing the same shape,
// but with only n points, uniformly distributed along the path from the start
// point to the end point.  Path should be in array-of-points format.
Sketchy.distributePointsAcrossPath = function (path, numberOfPoints) {
  const pathLength = Sketchy.computeLength(path, Sketchy.euclideanDistance);
  const delta = pathLength / numberOfPoints;
  let distanceCovered;
  let distanceToNextPoint;
  let angleToNextPoint;
  let nextPathIndex = 1;
  let currX = path[0].x;
  let currY = path[0].y;
  const result = [
    {
      x: currX,
      y: currY
    }
  ]; // Manually add the first point

  for (let i = 1; i < (numberOfPoints - 1); i++) {
    distanceCovered = 0;
    do {
      distanceToNextPoint = Sketchy.euclideanDistance(currX, currY, path[nextPathIndex].x, path[nextPathIndex].y);

      // Determine whether to jump to the next point or only move partially Last move
      // will occur in >= case (yes, it could happen in if or else)
      if (distanceToNextPoint <= delta - distanceCovered) {
        // Simply move to the next point
        currX = path[nextPathIndex].x;
        currY = path[nextPathIndex].y;
        nextPathIndex++;
        distanceCovered += distanceToNextPoint;
      } else {
        // Move partially
        angleToNextPoint = Math.atan2(path[nextPathIndex].y - currY, path[nextPathIndex].x - currX);
        currX += Math.cos(angleToNextPoint) * (delta - distanceCovered);
        currY += Math.sin(angleToNextPoint) * (delta - distanceCovered);
        distanceCovered = delta;
      }
    } while (distanceCovered < delta);
    // TODO: discretize currX and currY before pushing?
    result.push({ x: currX, y: currY });
  }
  // Manually add on the last point
  result.push(path[path.length - 1]);
  return result;
};

/* Betim's Algorithms */
// Compute the directed hausdorff distance of pixels1 and pixels2. Calculate the
// lowest upper bound over all points in shape1 of the distances to shape2.
// TODO: Make it faster!
Sketchy.hausdorff = function (points1, points2, vector2D) {
  let hMax = Number.MIN_VALUE;
  let hMin;
  let dist;
  for (let i = 0; i < points1.length; i++) {
    hMin = Number.MAX_VALUE;
    for (let j = 0; j < points2.length; j++) {
      dist = Sketchy.euclideanDistance(
        points1[i].x, points1[i].y, points2[j].x + vector2D.x, points2[j].y + vector2D.y
      );
      if (dist < hMin) {
        hMin = dist;
      } else if (dist === 0) {
        break;
      }
    }
    if (hMin > hMax) {
      hMax = hMin;
    }
  }
  return hMax;
};

// Compute hausdorffDistance hausdorff(shape1, shape2) and hausdorff(shape2,
// shape1) and return the maximum value.
Sketchy.hausdorffDistance = function (shape1, shape2, center1, center2) {
  const points1 = [];
  const points2 = [];
  const c1 = document.getElementById(shape1);
  const c2 = document.getElementById(shape2);
  const ctx1 = c1.getContext('2d');
  const ctx2 = c2.getContext('2d');
  const idata1 = ctx1.getImageData(0, 0, c1.width, c1.height);
  const idata2 = ctx2.getImageData(0, 0, c2.width, c2.height);
  for (let y1 = 0; y1 < c1.height; y1 += 4) {
    for (let x1 = 0; x1 < c1.width; x1 += 4) {
      if (idata1.data[(((x1 + y1) * c1.width) * 4) + 3] > 0) {
        points1.push({ x: x1, y: y1 });
      }
      if (idata2.data[(((x1 + y1) * c1.width) * 4) + 3] > 0) {
        points2.push({ x: x1, y: y1 });
      }
    }
  }
  const vector2D = {
    x: center1.x - center2.x,
    y: center1.y - center2.y
  };
  const h1 = Sketchy.hausdorff(points1, points2, vector2D);
  vector2D.x *= -1;
  vector2D.y *= -1;
  const h2 = Sketchy.hausdorff(points2, points1, vector2D);
  const accuracy = Math.max(h1, h2);
  return 1 - (((accuracy * Math.sqrt(2)) / 300) ** (1 / 1.4));
};

Sketchy.secondMoment = function (shape) {
  const c = document.getElementById(shape);
  const ctx = c.getContext('2d');
  const idata = ctx.getImageData(0, 0, c.width, c.height);
  const d = idata.data;
  let x;
  let y;
  const moment = {
    x: 0,
    y: 0
  };
  let tempX = 0;
  let tempY = 0;
  let size = 0;

  for (y = 0; y < c.height; y++) {
    for (x = 0; x < c.width; x++) {
      const value = d[(((x + y) * c.width) * 4) + 3];
      if (value > 0) {
        tempX += x;
        tempY += y;
        size++;
      }
    }
  }
  moment.x = tempX / size;
  moment.y = tempY / size;
  return moment;
};

export default Sketchy;
