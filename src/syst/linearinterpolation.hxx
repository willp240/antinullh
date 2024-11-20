//
// Taken from https://github.com/knotman90/linear_interpolation/blob/master/include/linearInterpolation.hxx
//

#ifndef SOLAR_ANALYSIS_LINEARINTERPOLATION_HXX
#define SOLAR_ANALYSIS_LINEARINTERPOLATION_HXX

#include <cassert>
#include <cstdio>

namespace FL {
    template<int DIM, typename T = long double>
    struct point {
        T coords[DIM];
        T val;

        inline T coord(const int c) const {
            assert(c >= 0  && c < DIM);
            return coords[c];
        }
    };


    template<class T>
    class LinearInterpolator {
    public:

        /*
             y
         ^         p1
         |        /
         |       /
         |      /
         |     p
         |    /
         |   p0
         |
         |
         o-------------------------> x
         */

        /*  P: point the lie between a and b
         *  a,b boundary of the 1D cuboid, a<b
         */
        static point<1, T>& Linear(point<1, T>& p,
                            const point<1, T>& a,
                            const point<1, T>& b,
                            int c = 0) {
            T x_d = (p.coord(c) - a.coord(c)) / (b.coord(c) - a.coord(c));

            p.val = Linear(a.val, b.val, x_d);
            return p;
        }

        /*  P: point the lie inside the cuboid defined by the first two values of v
         *     v[0] <= v[1]
         */
        static point<1, T>& Linear(point<1, T>& p, const point<1, T> *v, int c = 0) {
            T x_d = (p.coord(c) - v[0].coord(c)) / (v[1].coord(c) - v[0].coord(c));

            p.val = Linear(v[0].val, v[1].val, x_d);
            return p;
        }

        /*----------------BILINEAR------------------------
                   y
         ^
         |   p2.......p3
         |   .        .
         |   .        .
         |   .        .
         |   .        .
         |   p0.......p1
         |
         |
         o-------------------------> x

           p point that lie in the cuboid defined by the 4 values array v
           -------------------------------------------------*/

        static point<2, T>& Bilinear(point<2, T>& p, const point<2, T> *v) {
            T x_d = (p.coord(0) - v[0].coord(0)) / (v[1].coord(0) - v[0].coord(0));
            T y_d = (p.coord(1) - v[0].coord(1)) / (v[2].coord(1) - v[0].coord(1));

            p.val = Bilinear(v[0].val, v[1].val, v[2].val, v[3].val, x_d, y_d);

            return p;
        }
        /*----------------TRILINEAR------------------------
         *
                      y            p6____________p7
                      ^           / |            /|
                      |          /  |           / |
                      |         /   |          /  |
                      |        p2___|_________p3  |
                      |        |    |         |   |
                      |        |   p4_________|__p5
                      |        |   /          |  /
                      |        |  /           | /
                      |        | /            |/
                      |        p0_____________p1
                      |
                      |-------------------------> x
                      /
                     /
                    /
                   /
                  /
                z/


                p point that lies in the cuboid defined by the 8 values array v
        -------------------------------------------------*/

        static point<3, T>& Trilinear(point<3, T>& p, const point<3, T> *v) {
            T x_d = (p.coord(0) - v[0].coord(0)) / (v[1].coord(0) - v[0].coord(0));
            T y_d = (p.coord(1) - v[0].coord(1)) / (v[2].coord(1) - v[0].coord(1));
            T z_d = (p.coord(2) - v[0].coord(2)) / (v[4].coord(2) - v[0].coord(2));

            p.val = Trilinear(v[0].val,
                              v[1].val,
                              v[2].val,
                              v[3].val,
                              v[4].val,
                              v[5].val,
                              v[6].val,
                              v[7].val,
                              x_d,
                              y_d,
                              z_d);

            return p;
        }

    private:

        static inline T Linear(const T f0, const T f1, const T xd) {
            return f0 * (1.0 - xd) + f1 * xd;
        }

        static inline T Bilinear(const T f00,
                          const T f10,
                          const T f01,
                          const T f11,
                          const T xd,
                          const T yd) {
            const T c0 = f00 * (static_cast<T>(1.0) - xd) + f10 * xd;
            const T c1 = f01 * (static_cast<T>(1.0) - xd) + f11 * xd;

            return Linear(c0, c1, yd);
        }

        static inline T Trilinear(const T f000,
                           const T f100,
                           const T f010,
                           const T f110,
                           const T f001,
                           const T f101,
                           const T f011,
                           const T f111,
                           const T xd,
                           const T yd,
                           const T zd) {
            const T c0 = Bilinear(f000, f100, f010, f110, xd, yd);
            const T c1 = Bilinear(f001, f101, f011, f111, xd, yd);

            return Linear(c0, c1, zd);
        }
    }; // class interpolator
}

#endif //SOLAR_ANALYSIS_LINEARINTERPOLATION_HXX
