/**
 * Implementation of parametrized PSLQ algorithm in a complex case.
 * References: Analysis Of PSLQ, An Integer Relation Finding Algorithm,
 * Helaman Ferguson, David H. Bailey (1998)
 *
 * @param(tau) must lay between 1 and sqrt(2).
 * @example
 * var x = List.turnIntoCSList([ CSNumber.real(2), CSNumber.real(4) ]);
 * var pslq = new PSLQ(tau);
 * pslq.find(x) // Output will be a List too.
 */
var PSLQ = function(tau) {
    // rho = sqrt(2) for complex case
    // rho = 2 for real number field
    var rho = Math.sqrt(2);

    this.rho = CSNumber.real(rho);
    this.tau = tau && tau.ctype === "number" ? tau : CSNumber.real(Math.sqrt(4/3));

    if (this.tau.value.real <= 1 || this.tau.value.real >= rho) {
        throw new Error('tau must be bigger than 1 and smaller than sqrt(2)');
    }

    // tau, rho and gamma must satisfy: 1 / tau^2 = 1 / rho ^ 2 + 1 / gamma ^ 2
    this.gamma = CSNumber.inv(CSNumber.sqrt(CSNumber.real(1 / (this.tau.value.real * this.tau.value.real) - 1 / (rho * rho))));
};

/**
 * @param{x} input vector
 */
PSLQ.prototype.find = function(x) {
    var i, j;
    var result = [];
    var index;

    this.initialize(x);

    // Check whether x has zero component or not. If it has then relation is trivial.
    for (i = 0; i < this.n; i++) {
        if (CSNumber._helper.isAlmostZero(x.value[i])) {
            for (j = 0; j < this.n; j++) {
                result[j] = CSNumber.zero;
            }
            result[i] = CSNumber.real(1);

            return List.turnIntoCSList(result);
        }
    }

    while (index === void 0) {
        this.step1();
        this.step2();
        this.step3();
        index = this.shouldTerminate();
    }

    if (index === -1) {
        return false;
    }

    // Result is appearing as (n - 1) column of B (Lemma 5)
    for (j = 0; j < this.n; j++) {
        result.push(this.B.value[j].value[index]);
    }

    return List.turnIntoCSList(result);
};

/**
 * Initial step of PSLQ algorithm
 */
PSLQ.prototype.initialize = function(x) {
    // According to the definition 2, matrix Hx was introduced for vector x, |x| = 1.
    this.x = List.normalizeAbs(x);

    this.r = void 0;
    this.alpha = void 0;
    this.beta = void 0;
    this.lambda = void 0;
    this.delta = void 0;
    
    this.n = x.value.length;

    // Initialize A and B as identity matrices.
    this.A = List.idMatrix(CSNumber.real(this.n));
    this.B = List.idMatrix(CSNumber.real(this.n));

    this.computeHx();

    // attention: r is not passed
    this.hermiteReduction();
};

PSLQ.prototype.step1 = function() {
    var swap, R;

    this.choose_r();
    console.log(this.r);

    // Following parameters are used only in case r < n - 2 on step 2
    if (this.r < this.n - 2) {
        this.alpha = this.Hx.value[this.r].value[this.r];
        this.beta = this.Hx.value[this.r + 1].value[this.r];
        this.lambda = this.Hx.value[this.r + 1].value[this.r + 1];

        // delta = sqrt(beta * (beta*) + lambda * (lambda*));
        this.delta = CSNumber.sqrt(CSNumber.add(
            CSNumber.mult(this.beta, CSNumber.conjugate(this.beta)),
            CSNumber.mult(this.lambda, CSNumber.conjugate(this.lambda))));
    } 
    
    // Permutation matrix
    R = List.idMatrix(CSNumber.real(this.n));
    swap = R.value[this.r];
    R.value[this.r] = R.value[this.r + 1];
    R.value[this.r + 1] = swap;

    // Update x, Hx, A, B
    this.x = List.productVM(this.x, R);
    this.Hx = List.productMM(R, this.Hx);
    this.A = List.productMM(R, this.A);
    this.B = List.productMM(this.B, R);
};

PSLQ.prototype.choose_r = function() {
    var i, r;
    var n = this.n;
    var flag, cr, ci;
    var hAbs = [];

    for (i = 0; i < n - 1; i++) {
        hAbs.push(CSNumber.abs(this.Hx.value[i].value[i]));
    }

    // Find first r that:
    //     gamma^r * |h_rr| >= gamma^i * |h_ii|
    //     for all i: 0 <= i <= n - 2
    for (r = 0; r < n - 1; r++) {
        flag = true;

        cr = CSNumber.abs2(CSNumber.mult(CSNumber.powRealExponent(this.gamma, r), hAbs[r]));

        for (i = 0; i < n - 1; i++) {
            // TODO lazy calculation of power?
            ci = CSNumber.abs2(CSNumber.mult(CSNumber.powRealExponent(this.gamma, i), hAbs[i]));

            if (cr.value.real < ci.value.real) {
                flag = false;
            }
        }

        if (flag) {
            this.r = r;
            return;
        }
    }

    throw new Error('"r" was not found');
};

PSLQ.prototype.step2 = function() {
    if (this.r === this.n - 2) {
        // H remains unchanged
        return;
    }

    var i, m = this.r;

    // t0 = sqrt(H_mm^2 + H_(m+1)(m+1)^2)
    var t0 = CSNumber.sqrt(
        CSNumber.add(CSNumber.mult(this.Hx.value[m].value[m], this.Hx.value[m].value[m]),
            CSNumber.mult(this.Hx.value[m].value[m+1], this.Hx.value[m].value[m + 1])));

    // t1 = H_mm / t0
    var t1 = CSNumber.div(this.Hx.value[m].value[m], t0);

    // t2 = H_m(m+1) / t0
    var t2 = CSNumber.div(this.Hx.value[m].value[m + 1], t0);

    var t3, t4;

    for (i = m; i < this.n; i++) {
        // t3 = H_im
        t3 = this.Hx.value[i].value[m];
        // t4 = H_i(m+1)
        t4 = this.Hx.value[i].value[m + 1];

        // H_im = t1*t3 + t2*t4
        this.Hx.value[i].value[m] = CSNumber.add(CSNumber.mult(t1, t3), CSNumber.mult(t2, t4));
        // H_i(m+1) = t1*t4 - t2*t3
        this.Hx.value[i].value[m + 1] = CSNumber.sub(CSNumber.mult(t1, t4), CSNumber.mult(t2, t3));
    }
   // List.println(this.Hx);
};

PSLQ.prototype.step3 = function() {
    this.hermiteReduction();
};

/**
 * Checks whether algorithm should be terminated.
 * Returns true if one of the foolowing conditions is held:
 * - there is zero component in x vector;
 * - there is zero component on the diagonal of Hx
 * - any entry of A matrix exceeds the level of numeric precision being used
 */
PSLQ.prototype.shouldTerminate = function() {
    var i, j;

    for (i = 0; i < this.n - 1; i++) {
        if (CSNumber._helper.isAlmostZero(this.x.value[i])) {
            return i;
        }

        if (CSNumber._helper.isAlmostZero(this.Hx.value[i].value[i])) {
            return -1;
        }
    }

    if (CSNumber._helper.isAlmostZero(this.x.value[this.n - 1])) {
        return this.n - 1;
    }

    return void 0;
};

/**
 * Computes Hx, n x (n - 1) matrix.
 */
PSLQ.prototype.computeHx = function() {
    var s = [];
    var i, j, k;
    var Hx = [];
    var xj2 = [];

    for (j = 0; j < this.n; j++) {
        xj2.push(CSNumber.mult(this.x.value[j], CSNumber.conjugate(this.x.value[j])));
    }

    // Define the partial sums s_j
    for (j = 0; j < this.n; j++) {
        s[j] = CSNumber.zero;

        for (k = j; k < this.n; k++) {
            s[j] = CSNumber.add(s[j], xj2[k]);
        }

        s[j] = CSNumber.sqrt(s[j]);
    }

    for (i = 0; i < this.n; i++) {
        // add empty row
        Hx[i] = [];

        for (j = 0; j < i; j++) {
            // Compute Hx using the following formula: Hx_ij = -(x_i*) * x_j / (s_j * s_(j + 1))
            Hx[i][j] = CSNumber.sub(CSNumber.zero, CSNumber.div(
                CSNumber.mult(CSNumber.conjugate(this.x.value[i]), this.x.value[j]), CSNumber.mult(s[j], s[j + 1])));
        }

        if (i < this.n - 1) {
            // Hx_ii = s_(i+1) / s_i
            Hx[i][i] = CSNumber.div(s[i + 1], s[i]);


            for (j = i + 1; j < this.n - 1; j++) {
                Hx[i][j] = CSNumber.zero;
            }
        }

        Hx[i] = List.turnIntoCSList(Hx[i]);
    }

    this.Hx = List.turnIntoCSList(Hx);
};

/**
 * Computes modified Hermite reduction (from section 8: Computer implementation).
 * Updates matrices Hx, A, B and vector x.
 */
PSLQ.prototype.hermiteReduction = function() {
    var i, j, k;
    var isInitStep = this.r === void 0;
    var t, jmax;
    var i0 = isInitStep ? 1 : this.r + 1;
//    List.println(this.Hx);

    for (i = i0; i < this.n; i++) {
        jmax = isInitStep ? i - 1 : Math.min(i - 1, this.r + 1);

        for (j = jmax; j >= 0; j--) {
            t = CSNumber.round(CSNumber.div(this.Hx.value[i].value[j], this.Hx.value[j].value[j]));
            // Compute formula: x_j = x_j - t * x_i
            this.x.value[j] = CSNumber.add(this.x.value[j], CSNumber.mult(t, this.x.value[i]));

            for (k = 0; k <= j; k++) {
                // H_ik = H_ik - t * H_jk
               // console.log('t = ' + t.value.real, ' i = ' + i, ' j = ' + j, ' k = ' + k);
                this.Hx.value[i].value[k] = CSNumber.sub(this.Hx.value[i].value[k], CSNumber.mult(t, this.Hx.value[j].value[k]));
               // List.println(this.Hx.value[2]);
            }

            for (k = 0; k < this.n; k++) {
                // A_ik = A_ik - t * A_jk
                this.A.value[i].value[k] = CSNumber.sub(this.A.value[i].value[k], CSNumber.mult(t, this.A.value[j].value[k]));

                // B_kj = B_kj + t * B_ki
                this.B.value[k].value[j] = CSNumber.add(this.B.value[k].value[j], CSNumber.mult(t, this.B.value[k].value[i]));
            }
        }
    }
};
