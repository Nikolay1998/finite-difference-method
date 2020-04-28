public class Solver {

  private Double paramC = 1.65;
  private Double paramK = 0.59;
  private Double T = 30.0;
  private Double R = 4.0;
  private Integer K = 1000;
  private Integer I = 1000;

  Solver() {
    Double[][] u = new Double[K][I];

    Double[] b = new Double[I - 2];
    Double[] a = new Double[I - 2];
    Double[] c = new Double[I - 2];
    Double[] d = new Double[I - 2];


    Double hr = R / I;
    Double ht = T / K;

    u[0] = zeroLayerCalc(hr);



    Double coef = ht / paramC;
    Integer matrSize = I - 2;

    Double r = 0.0;
    r += hr;
    for (int i = 0; i < matrSize; i++) {
      a[i] = (paramK / (r * hr) - paramK / (hr * hr)) * coef;
      r += hr;
    }


    r = 0.0;
    r += hr;
    for (int i = 0; i < matrSize; i++) {
      c[i] = (-paramK / (r * hr) - paramK / (hr * hr)) * coef;
      r += hr;
    }


    r = hr;
    Double b1 = 1.0 + (2 * paramK * ht) / (paramC * hr * hr);
    b[0] = b1 + a[0] / ((paramC * hr * hr / (ht * 6 * paramK)) + 1);

    for (int i = 1; i < matrSize; i++) {
      b[i] = b1;
    }

    b[matrSize-1] = b1 + c[matrSize-1] / ((hr * hr * paramC / (ht * 2 * paramK)) + 1);

    Double[] shuttleResult;

    for (int k = 1; k < K - 1; k++) {
      d[0] = u[k - 1][1] - a[0] * u[k - 1][0] / (1 + (6 * paramK * ht / (hr * hr * paramC)));
      for (int i = 1; i < matrSize - 1; i++) {
        d[i] = u[k - 1][i + 1];
      }
      //System.out.println(u[k-1][I-1-1]);
      d[matrSize-1] = u[k-1][I-1-1] - c[matrSize-1]/(1+ (2*paramK*ht/(paramC * hr * hr)));
      shuttleResult = shuttle(b, a, c, d);
      u[k][0] = (shuttleResult[0] / ((paramC * hr * hr / (ht * 6 * paramK)) + 1)) + (u[k - 1][0] / (1 + (6 * paramK * ht / (hr * hr * paramC))));
      u[k][I-1] = shuttleResult[matrSize-1]/(1 + (hr*hr*paramC/(ht*2*paramK))) + (u[k-1][I-1] / (1 + 2*paramK*ht/(paramC*hr*hr)));
      for (int i = 0; i < matrSize; i++) {
        u[k][i + 1] = shuttleResult[i];
      }
    }

    for (int i = 0; i < I; i++) {
      System.out.println(i + ") " + u[998][i]);
    }

  }

  private Double[] zeroLayerCalc(Double hr) {
    Double r = 0.0;
    Double[] result = new Double[I];
    for (int i = 0; i < I; i++) {
      result[i] = psi(r);
      r += hr;
    }
    return result;
  }

  private Double[] shuttle(Double[] b, Double[] a, Double[] c, Double d[]) {
    Integer n = d.length;
    Double[] y = new Double[n];
    Double[] alpha = new Double[n];
    Double[] betta = new Double[n];
    Double[] x = new Double[n];
    y[0] = b[0];
    alpha[0] = -c[0] / y[0];
    betta[0] = d[0] / y[0];
    for (int i = 1; i < n - 1; i++) {
      y[i] = b[i] + a[i] * alpha[i - 1];
      alpha[i] = -c[i] / y[i];
      betta[i] = (d[i] - a[i] * betta[i - 1]) / y[i];
    }
    y[n - 1] = b[n - 1] + a[n - 1] * alpha[n - 1 - 1];
    betta[n - 1] = (d[n - 1] - a[n - 1] * betta[n - 1 - 1]) / y[n - 1];

    x[n - 1] = betta[n - 1];
    for (int i = n - 1 - 1; i >= 0; i--) {
      x[i] = alpha[i] * x[i + 1] + betta[i];
    }

    return x;
  }


  private Double psi(Double r) {
    return Math.cos(Math.PI * r / (2 * R));
  }
}
