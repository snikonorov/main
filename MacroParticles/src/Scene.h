#ifndef SCENE_H_INCLUDED
#define SCENE_H_INCLUDED

#include "SGL.h"

/// basic force magnitude ('regularized' Lennard-Jones)
///
decltype(auto) A0(double R0, double d, double g = 0.5)
{
    const double C0 = 8e3;

    ///const double g = 0.5;
    ///const double g = 0.9;       /// good for creating holographic patterns
    ///const double g = 0.3;

    const double d1 = g*R0;
    const double d2 = g*R0;

    if(d <= R0 - d1)
    {
        return -C0;
    }
    else if(d > R0 - d1 && d <= R0)
    {
        double q = (d - R0)/d1;

        return C0*(R0 - d)*(q*q - 3)/(2*d1);
    }
    else if(d > R0 && d <= R0 + d2)
    {
        double q = (R0 + d2 - d)/d2;

        return 3*C0*(d - R0)*q*q/(2*d1);
    }
    else
    {
        return 0.0;
    }
}

/// basic spring force
/*decltype(auto) A0(double R0, double d)
{
    double a = 3;
    return d > a*R0 ? 0.0 : (d < R0 ? 2e4 : 2e2)*(d/R0 - 1);
}*/

/// order 5 asymmetric polynomial
/*decltype(auto) A0(double R0, double d, double a = 1)
{
    const double C0 = 1e2;

    a = 1.0;

    double r = d/R0;

    if(r <= 4)
    {
        return -C0*(32*(a - r)*(r - 2)*(r - 2)*r*r)/(18*a - 9);
    }
    else
    {
        return 0.0;
    }
}*/

/// rational approximation of a simplified regularized LG
///
/*decltype(auto) A0(double R0, double d, double a = 1)
{
    const double C0 = 8e3;

    a = 1.0;

    double r = d/R0;

    return C0*(r - 1)/(1 - r*(1 - r*(1 - r*(1 - r*(1 - r*(1 - 10*r))))));
}*/

///----------------------------------------------------------------

struct Object
{
    virtual double operator()(Point r) const
    {
        return 0.0;
    }

    virtual void render(Image& canvas) const
    {
    }

    virtual void update(double dt)
    {
    }
};

/// volume particle
struct VParticle
{
    Point p;
    Point v;
    Point a;
    double R0;

    double g;       /// 'adhesive' index

    pixel color;

    VParticle(Point p, double R0 = 20, double g = 0.5, const pixel& color = pixel(0, 0, 220))
             : p(p), v(0.0), a(0.0), R0(R0), g(g), color(color)
    {
    }

    /// particle `density distribution` function
    ///
    double operator()(Point r) const
    {
        //double q = (p - r).Sqr()/(R0*R0);
        ///double q = 0.11*(p - r).Sqr()/(R0*R0);

        //return exp(-q*q);
        ///return q > 1 ? 0.0 : (1 - q)*(1 - q);

        const double d = 0.1;               /// fall-off range length

        double a = (p - r).Norm()/(1.2*R0) - 1;
        double q = a/d;
        double q2 = q*q;

        return a > 0 ? 0.0 : -q*q2*(10 + 15*q + 6*q2);
    }
};

///----------------------------------------------------------------

/// particle system
struct PSystem: public Object, public std::vector<VParticle>
{
    PSystem() = default;

    PSystem(const std::vector<VParticle>& data) : std::vector<VParticle>(data)
    {
    }

    double operator()(Point r) const override
    {
        //const double d_max = 0.5;

        double d = sum(map([r](const auto& P){ return P(r); }, *this));

        return d >= 1 ? 1.0 : d;
    }

    void update(double dt = 1e-2) override
    {
        const double G0 = 10e0; /// ~friction

        ///------------------------------------------------

        unsigned const N = 10;
        dt /= N;

        for(unsigned n = 0; n < N; n++)
        {
            for(unsigned k = 0; k < size(); k++)
            {
                auto& P1 = at(k);

                P1.a = Point(0.0);

                for(unsigned p = 0; p < size(); p++)
                {
                    if(k != p)
                    {
                        auto& P2 = at(p);

                        Point dp = P2.p - P1.p;
                        double d = dp.Norm();

                        double R0 = (P2.R0 + P1.R0);

                        if(d > 1)
                        {
                            P1.a += A0(R0, d, 0.5*(P1.g + P2.g))*dp/d;
                        }
                    }
                }
            }

            for(unsigned k = 0; k < size(); k++)
            {
                auto& P = at(k);

                ///P.v += (P.a - G0*P.v*P.v.Norm())*dt;
                P.v += (P.a - G0*P.v)*dt;

                P.p += P.v*dt;
            }
        }
    }

    void render(Image& canvas) const override
    {
        /// draw particle 'density'
		///
        if(0)
        {
            for(int x = 0; x < canvas.Width; x++)
            {
                for(int y = 0; y < canvas.Height; y++)
                {
                    double I = (*this)(Point(x, y));

                    canvas.PutPixel(x, y, pixel(255*I));
                }
            }
        }

        ///----------------------------------------------------

        /// draw connection graph
        if(1)
        {
            for(unsigned k = 0; k < size(); k++)
            {
                for(unsigned p = k+1; p < size(); p++)
                {
                    double dp = (at(k).p - at(p).p).Norm();
                    double R0 = at(k).R0 + at(p).R0;

                    if(dp <= R0 + 0.5)
                    {
                        /// ~ strain
                        double S = std::max(0.0, 1 - dp/R0);

                        double a = ceil(255*std::min(1.0, 4*S));

                        /// 'strain line'
                        DrawDDALine
                        (
                            at(k).p, at(p).p,
                            pixel(a, 255, 0, 0)
                        );
                    }
                }
            }
        }

        ///----------------------------------------------------

        for(const auto& P: *this)
        {
            double a = 1;
            const pixel& C = P.color;

            canvas.FillCircle(P.p.x, P.p.y, 3, 0.4*a*C);
            canvas.DrawDDACircle(P.p.x, P.p.y, 3, a*C);

            //canvas.DrawDDACircle(P.p.x, P.p.y, P.R0, pixel(255, 0, 0));
            canvas.DrawDDACircle(P.p.x, P.p.y, P.R0, 0.2*a*C);
        }
    }
};

struct Scene: public TComponent
{
    Image canvas;
    std::vector<Object*> data;

    ///-------------------------------------------

    void render(double dt = 0)
    {
        canvas.resize(AppSize.x, AppSize.y);
        canvas.Clear();

        for(auto* obj: data)
        {
            obj->update(dt);
            obj->render(canvas);
        }
    }

    void ProcessEvent(sEvent& event)
    {
        static clock_t T = 0;
        static bool freeze = false;
        static sFont* font = new sFont(14);

        using namespace sEventCodes;

        switch(event)
        {
            case Keyboard_Down:

                switch(sKey.Key)
                {
                    case 'F':
                        freeze = !freeze;
                        break;
                }
                break;

            case Window_Paint:
                {
                    double dt = 0;

                    if(T)
                        dt = 1.*(clock() - T)/CLOCKS_PER_SEC;

                    T = clock();

                    render(freeze ? 0.0 : dt);
                    DrawImage(canvas, 0, 0);

                    if(freeze)
                    {
                        Text_out(font, 10, 10, "Frozen", pixel(100, 0, 0, 255));
                    }
                }
                break;

            default:
                break;
        }
    }
};

#endif // SCENE_H_INCLUDED
