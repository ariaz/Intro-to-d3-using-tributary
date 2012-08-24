/**
 * Tuebingen demo.
 *
 * Based on Greg Egan's applet demonstrating projections of a 5-dimensional
 * lattice onto a 2D plane.
 */
var Tuebingen = function() {
  this.bfac0 = 0.2;
}

var random = function() {
  return Math.random();
}

var mainParameter = -1;

var dim = function() {
  var len = arguments.length;
  if (len > 0) {
    var n = arguments[0];
    var a = new Array(n);
    if (len > 1) {
      var args = [];
      for (var i=1; i<len; i++) {
        args.push(arguments[i]);
      }
      for (var i=0; i<n; i++) {
        a[i] = dim.call(this, args);
      }
    } else {
      for (var i=0; i<n; i++) {
        a[i] = 0;
      }
    }
    return a;
  }
}

Tuebingen.prototype.setup = function(l, i) {
  var minD = 4, maxD = 14, numD = (maxD - minD) + 1;
  if (i > 0) mainParameter = i;
  else if (mainParameter < 0) for (; mainParameter == 5 || mainParameter < 0; mainParameter = minD + Math.floor(numD * random()));
  var dA = mainParameter;
  var n1 = minD + ((dA - minD) + 1) % numD;
  var n2 = minD + (((dA - minD) + numD) - 1) % numD;
  var dI = this.DI = dA + 1;
  var oddDim = this.oddDim = dI % 2 != 0;
  var bx = this.bx = dim(dI);
  var by = this.by = dim(dI);
  var p = this.p = dim(dI);
  var phase = 6 * random();
  var sA = Math.PI / dI, hA = sA / 2;
  var theta = 2 * sA;
  var norm = Math.sqrt(2 / dI);
  var sp = 0;
  for (var k = 0; k < dI; k++) {
    bx[k] = norm * Math.cos(phase + theta * k);
    by[k] = norm * Math.sin(phase + theta * k);
    sp += p[k] = random();
  }

  sp /= dI;
  for (var k = 0; k < dI; k++)
    p[k] -= sp;

  if (oddDim)
    marg = (0.5 * norm) / Math.sin(hA);
  else
    marg = norm / Math.sin(sA);
  var dpSpan = this.dpSpan = dim(dA);
  var dx = this.dx = dim(dA);
  var dy = this.dy = dim(dA);
  var dp = this.dp = dim(dA);
  for (var i = 0; i < dA; i++) {
    var j2 = dA - i;
    dpSpan[i] = dim(j2, 2);
    var ad = dx[i] = dim(j2);
    var ad1 = dy[i] = dim(j2);
    var ad2 = dp[i] = dim(j2);
    for (var j = i+1, k = 0; j < dI; j++, k++) {
      ad[k] = this.bx[i] - this.bx[j];
      ad1[k] = by[i] - by[j];
      ad2[k] = p[i] - p[j];
    }
  }

  var ntrip = this.ntrip = Math.floor((dI * dA * (dA - 1)) / 6);
  tripCoords = dim(ntrip, dI);
  mapCoords = dim(ntrip, dI);
  var solAx = this.solAx = dim(ntrip);
  var solBx = this.solBx = dim(ntrip);
  var solCx = this.solCx = dim(ntrip);
  var solAy = this.solAy = dim(ntrip);
  var solBy = this.solBy = dim(ntrip);
  var solCy = this.solCy = dim(ntrip);
  var l = 0;
  for (var i = 0; i < dI - 2; i++) {
    for (var j = i + 1; j < dA; j++) {
      for (var k = j + 1; k < dI; k++) {
        var tc = tripCoords[l], mc = mapCoords[l];
        tc[0] = i;
        tc[1] = j;
        tc[2] = k;
        mc[i] = 0;
        mc[j] = 1;
        mc[k] = 2;
        var m = 3;
        for (var j5 = 0; j5 < dI; j5++) {
          if (j5 != i && j5 != j && j5 != k) {
            tc[m] = j5;
            mc[j5] = m;
            m++;
          }
        }
        var denom = bx[k] * (by[i] - by[j]) + bx[i] * (by[j] - by[k]) + bx[j] * (by[k] - by[i]);
        solAx[l] = (by[i] - by[k]) / denom;
        solBx[l] = (by[j] - by[i]) / denom;
        solCx[l] = (by[k] * (p[i] - p[j]) + by[i] * (p[j] - p[k]) + by[j] * (p[k] - p[i])) / denom;
        solAy[l] = (bx[k] - bx[i]) / denom;
        solBy[l] = (bx[i] - bx[j]) / denom;
        solCy[l] = (bx[k] * (p[j] - p[i]) + bx[i] * (p[k] - p[j]) + bx[j] * (p[i] - p[k])) / denom;
        l++;
      }
    }
  }

  // Determine tile colours
  var Tshape = this.Tshape = dim(ntrip);
  var Torient = this.Torient = dim(ntrip);
  var cfac = this.cfac = random();
  var maxDet = 0;
  for (var tt = 0; tt < ntrip; tt++) {
    var tc = tripCoords[tt];
    var i = tc[0];
    var j = tc[1];
    var k = tc[2];
    var ij = j - i - 1, ik = k - i - 1;
    var xe = dx[i][ij], ye = dy[i][ij];
    var xf = dx[i][ik], yf = dy[i][ik];
    var det = Math.sqrt(Math.abs(xf * ye - yf * xe));
    Tshape[tt] = det;
    if (det > maxDet)
      maxDet = det;
    Torient[tt] = Math.abs((Math.atan2(xf, yf) + Math.atan2(xe, ye)) / (2 * Math.PI));
  }

  for (var j3 = 0; j3 < ntrip; j3++)
    Tshape[j3] /= maxDet;
}

Tuebingen.prototype.draw = function(i, j, k, yoffs, flag) {
  var x, y;
  var cfac = this.cfac;
  var lines = [];
  var line = [];
  var dx = this.dx, dy = this.dy, dp = this.dp, dpSpan = this.dpSpan,
      p = this.p, bx = this.bx, by = this.by,
      ntrip = this.ntrip,
      solAx = this.solAx, solBx = this.solBx, solCx = this.solCx,
      solAy = this.solAy, solBy = this.solBy, solCy = this.solCy,
      Tshape = this.Tshape, Torient = this.Torient;
  var wc = i / 2;
  var h2 = (k - j) / 2;
  var hc = j + h2;
  if (flag) {
    thisSeed = random();
    this.setup(thisSeed, -1);
  }
  var dA = mainParameter;
  var dI = this.DI;
  var uuu = dim(dA);
  var vvv = dim(3, dI);
  var rp = dim(dI);
  var span0 = 10 / Math.sqrt(dI / 5);
  var inc = span0 / 300;
  var span = inc * i;
  var xm1 = wc * inc + marg, xm0 = -xm1;
  var ym1 = (yoffs + Math.floor(h2)) * inc + marg;
  var ym0 = (yoffs - Math.floor(h2)) * inc - marg;
  var xm = [xm0, xm1];
  var ym = [ym0, ym1];
  var xep = dim(2);
  var yep = dim(2);
  for (var i = 0; i < dA; i++) {
    var dxi = dx[i], dyi = dy[i], dpi = dp[i];
    for (var j = i + 1, k = 0; j < dI; j++, k++) {
      var lc = Number.MAX_VALUE;
      var uc = Number.MIN_VALUE;
      for (var jj = 0; jj < 2; jj++) {
        for (var ii = 0; ii < 2; ii++) {
          var dot = Math.ceil(dxi[k] * xm[ii] + dyi[k] * ym[jj] + dpi[k]);
          if (dot < lc) lc = dot;
          if (--dot > uc) uc = dot;
        }
      }
      dpSpan[i][k][0] = lc;
      dpSpan[i][k][1] = uc;
    }
  }

  for (var tt = 0; tt < ntrip; tt++) {
    var tc = tripCoords[tt];
    var mc = mapCoords[tt];
    var sax = solAx[tt], sbx = solBx[tt], scx = solCx[tt];
    var say = solAy[tt], sby = solBy[tt], scy = solCy[tt];
    var i = tc[0];
    var j = tc[1];
    var k = tc[2];
    var ij = j - i - 1;
    var ik = k - i - 1;
    var dxi = dx[i], dyi = dy[i]; var dpi = dp[i];
    var xij = dxi[ij], xik = dxi[ik];
    var yij = dyi[ij], yik = dyi[ik];
    var pij = dpi[ij], pik = dpi[ik];
    var tshape = Tshape[tt], torient = Torient[tt];
    for (var ijD = dpSpan[i][ij][0]; ijD <= dpSpan[i][ij][1]; ijD++) {
      var ae = ijD - pij;
      var nn = 0;
      for (var ii = 0; ii < 2 && nn < 2; ii++) {
        x = xm[ii];
        y = (ae - xij * x) / yij;
        if (y >= ym0 && y <= ym1) {
          xep[nn] = x;
          yep[nn++] = y;
          if (nn == 2)
            break;
        }
        y = ym[ii];
        x = (ae - yij * y) / xij;
        if (x >= xm0 && x <= xm1) {
          xep[nn] = x;
          yep[nn++] = y;
        }
      }

      var lc0 = Math.ceil(xik * xep[0] + yik * yep[0] + pik);
      var uc0 = Math.floor(xik * xep[1] + yik * yep[1] + pik);
      if (uc0 < lc0) {
        var tmp = lc0;
        lc0 = uc0 + 1;
        uc0 = tmp - 1;
      }
      var xp = sax * ijD + scx;
      var yp = say * ijD + scy;
      for (var ikD = lc0; ikD <= uc0; ikD++) {
        x = xp + sbx * ikD;
        y = yp + sby * ikD;
        for (var k7 = 0; k7 < dI; k7++)
          rp[k7] = p[k7] + x * bx[k7] + y * by[k7];

        var sum = ijD + ikD;
        uuu[0] = ijD;
        uuu[1] = ikD;
        var x0 = rp[i];
        for (var i8 = 2; i8 < dA; i8++)
          sum += uuu[i8] = Math.floor(x0 - rp[tc[i8 + 1]]);

        var k8 = sum <= 0 ? dI - -sum % dI : sum % dI;
        var neg = k8 == 1;
        if (!neg && k8 != 2)
          continue;
        var q;
        if (neg) {
          for (var l8 = 0; l8 < dA; l8++)
            uuu[l8] = -uuu[l8] - 1;

          q = -(sum - 1) / dI - 1;
        } else {
          uuu[0]--;
          uuu[1]--;
          q = (sum - 2) / dI;
        }
        for (var m = 0; m < dI; m++) {
          var mcm = mc[m] - 1;
          vvv[0][m] = vvv[1][m] = mcm >= 0 ? q - uuu[mcm] : q;
        }

        vvv[1][i]++;
        for (var j9 = 0; j9 < dI; j9++)
          vvv[2][j9] = vvv[1][j9];

        vvv[1][j]--;
        vvv[2][k]--;
        var sign = neg ? -1 : 1;
        var yy = 0;
        for (var m = 0; m < 3; m++) {
          x = 0; y = 0;
          for (var mm = 0; mm < dI; mm++) {
            var delta = (sign * vvv[m][mm]) - p[mm];
            x += bx[mm] * delta;
            y += by[mm] * delta;
          }
          var xc = x / inc;
          var yc = y / inc - yoffs;
          if (m == 0) {
            line = [];
            yy = y;
          }
          line.push([xc, yc]);
        }
        line.push(line[0]);
        var sfac1 = (1 - Math.sin(yy / 171)) / 2;
        var bfac = this.bfac0 + (1 - this.bfac0) * sfac1;
        line.c = d3.hsl(Math.ceil(360 * ((cfac + tshape + yy / 111) % 1)),
          0.7 * sfac1 + 0.3, bfac * torient + (1 - bfac)).rgb();
        lines.push(line);
      }
    }
  }
  return lines;
}
