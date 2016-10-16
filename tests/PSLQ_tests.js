var rewire = require("rewire");
var cindyJS = rewire('../build/js/exposed.js');

var List = cindyJS.__get__('List');
var CSNumber = cindyJS.__get__('CSNumber');
var PSLQ = cindyJS.__get__('PSLQ');
var expect = require("chai").expect;

var eps = 1e-15;

describe('PSLQ', function() {
    var pslq, tau, expected, result;

    beforeEach(function() {
        pslq = new PSLQ();
    });

    describe('initialization of PSLQ object', function() {
        afterEach(function() {
            CSNumber._helper.isAlmostEqual(result, expected, eps).should.be.true;
        });

        it('gamma', function() {
            result = pslq.gamma;
            expected = CSNumber.real(2);
        });

        it('tau', function() {
            result = pslq.tau;
            expected = CSNumber.real(Math.sqrt(4/3));
        });

        it('rho', function() {
            result = pslq.rho;
            expected = CSNumber.real(Math.sqrt(2));
        });
    });

    describe("initialize", function() {
        beforeEach(function() {
            pslq.initialize(List.realVector([1, 2, 4]));
        });

        describe("matrices and vectors", function() {
            afterEach(function() {
                List.almostequals(result, expected).value.should.be.true;
            });

            it("A", function() {
                expected = List.realMatrix([
                    [1, 0, 0],
                    [0, 1, 0],
                    [0, 0, 1]
                ]);
                result = pslq.A;
            });

            it("B", function() {
                expected = List.realMatrix([
                    [1, 0, 0],
                    [0, 1, 0],
                    [0, 0, 1]
                ]);
                result = pslq.B;
            });

            it("x", function() {
                var t = Math.sqrt(21);

                expected = List.realVector([1 / t, 2 / t, 4 / t]);
                result = pslq.x;
            });
        });

        describe("constants", function() {
            it("n", function() {
                pslq.n.should.be.equal(3);
            });

            it("r", function() {
                expect(pslq.r).to.be.undefined;
            });

            it("alpha", function() {
                expect(pslq.alpha).to.be.undefined;
            });

            it("beta", function() {
                expect(pslq.beta).to.be.undefined;
            });

            it("lambda", function() {
                expect(pslq.lambda).to.be.undefined;
            });

            it("delta", function() {
                expect(pslq.delta).to.be.undefined;
            });
        });
    });

    describe('Compute Hx', function() {
        afterEach(function() {
            List.almostequals(result, expected).value.should.be.true;
        });

        it('real x = [ 1 2 4 ]', function() {
            pslq.x = List.normalizeAbs(List.turnIntoCSList([CSNumber.real(1), CSNumber.real(2), CSNumber.real(4)]));
            pslq.n = 3;

            pslq.computeHx();
            result = pslq.Hx;
            expected = List.realMatrix([
                [ Math.sqrt(20/21), 0 ],
                [ -1 / Math.sqrt(105), 4 / Math.sqrt(20) ],
                [ -2 / Math.sqrt(105), - 1 / Math.sqrt(5) ]
            ]);
        });

        it('complex x = [ i 2i 4i ]', function() {
            pslq.x = List.normalizeAbs(List.turnIntoCSList([CSNumber.complex(0, 1), CSNumber.complex(0, 2), CSNumber.complex(0, 4)]));
            pslq.n = 3;

            pslq.computeHx();
            result = pslq.Hx;
            expected = List.realMatrix([
                [ Math.sqrt(20/21), 0 ],
                [ - 1 / Math.sqrt(105), 4 / Math.sqrt(20) ],
                [ - 2 / Math.sqrt(105), - 1 / Math.sqrt(5) ]
            ]);
        });
    });

    it("r", function() {
        pslq.initialize(List.turnIntoCSList([CSNumber.real(1), CSNumber.real(2), CSNumber.real(4)]));

        pslq.choose_r();
        pslq.r.should.be.equal(1);
    });

    describe("Hermite reduction", function() {
        beforeEach(function() {
            pslq.x = List.normalizeAbs(List.realVector([1, 2, 4]));

            pslq.r = void 0;
            pslq.n = 3;

            // Initialize A and B as identity matrices.
            pslq.A = List.idMatrix(CSNumber.real(3));
            pslq.B = List.idMatrix(CSNumber.real(3));
        });

        afterEach(function() {
            List.almostequals(result, expected, eps).value.should.be.true;
        });


        it('A', function() {
            result = pslq.A;
            expected = List.realMatrix([
                [ 1, 0, 0 ],
                [ 0, 1, 0 ],
                [ 0, 0, 1 ]
            ]);
        });

        it('B', function() {
            result = pslq.B;
            expected = List.realMatrix([
                [ 1, 0, 0 ],
                [ 0, 1, 0 ],
                [ 0, 0, 1 ]
            ]);
        });

        it('x', function() {
            result = pslq.x;
            expected = List.normalizeAbs(List.turnIntoCSList([CSNumber.real(1), CSNumber.real(2), CSNumber.real(4)]));
        });
    });

    describe("step 1", function() {
        describe("vectors and matrices", function() {
            beforeEach(function() {
                pslq.initialize(List.realVector([1, 2, 4]));
                pslq.step1();
            });

            afterEach(function() {
                List.almostequals(result, expected, eps).value.should.be.true;
            });

            it("A", function() {
                result = pslq.A;
                expected = List.realMatrix([
                    [1, 0, 0],
                    [0, 0, 1],
                    [0, 1, 0]
                ]);
            });

            it("B", function() {
                result = pslq.B;
                expected = List.realMatrix([
                    [1, 0, 0],
                    [0, 0, 1],
                    [0, 1, 0]
                ]);
            });

            it("H", function() {
                result = pslq.Hx;
                expected = List.realMatrix([
                    [ Math.sqrt(20/21), 0 ],
                    [ -2 / Math.sqrt(105), - 1 / Math.sqrt(5) ],
                    [ -1 / Math.sqrt(105), 4 / Math.sqrt(20) ]
                ]);
            });

            it("x", function() {
                result = pslq.x;
                var t = Math.sqrt(21);
                expected = List.realVector([1/t, 4/t, 2/t]);
            });
        });

   //     describe("constants", function() {
   //         beforeEach(function() {
   //             pslq.initialize(List.realVector([1, 2, 4]));
   //             pslq.step1();
   //             pslq.step2();
   //             pslq.step3();
   //             List.println(pslq.Hx);
   //             console.log(pslq.r);
   //             pslq.step1();
   //         });

   //         it("alpha", function() {
   //             expect(pslq.alpha).to.be.undefined;
   //         });

   //         it("beta", function() {
   //             expect(pslq.beta).to.be.undefined;
   //         });

   //         it("lambda", function() {
   //             expect(pslq.lambda).to.be.undefined;
   //         });
   //     });
    });

    describe("step 2", function() {
        beforeEach(function() {
            pslq.initialize(List.realVector([1, 2, 4]));
            pslq.step1();
            pslq.step2();
        });

        afterEach(function() {
            List.almostequals(result, expected, eps).value.should.be.true;
        });


    });
});
