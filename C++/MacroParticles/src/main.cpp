#include "Scene.h"

int main()
{
    auto* w = SCreateWindow(1000, 800, L"MacroParticles");

    auto* P = new PSystem();
    auto* s = new Scene();

    s->data.push_back(P);
    ///s->render();

    *w += s;

    /// particle control component
    ///
    *w += [P](sEvent& event)
    {
        static unsigned n = 0;
        const double R0 = 20;

        using namespace sEventCodes;

        switch(event)
        {
            case Mouse_LeftClick:
                if(sKey.Control)
                {
                    P->push_back(VParticle(sMouse.Pos(), R0));
                }
                break;

            case Mouse_Move:
                {
                    if(n && sKey.Shift)
                    {
                        P->at(n-1).p = sMouse.Pos();
                    }
                }
                break;

            case Keyboard_Down:

                switch(sKey.Key)
                {
                    case 8:
                        if(P->size())
                            P->pop_back();
                        break;

                    case 46:
                        P->clear();
                        break;
                }

                if(sKey.Shift)
                {
                    if(!n)
                    {
                        P->push_back(VParticle(sMouse.Pos()));
                        n = P->size();
                    }
                }

                break;

            case Keyboard_Up:
                {
                    if(!sKey.Shift)
                    {
                        if(n)
                        {
                            P->erase(P->begin() + (n - 1));
                            n = 0;
                        }
                    }
                }
                break;

            /*case Mouse_Wheel:
                for(auto& p: *P)
                {
                    p.R0 *= sMouse.Wheel > 0 ? 1.1 : 1/1.1;
                }

                s.render();
                break;*/

            default:
                break;
        }
    };

    ///------------------------------------------------------------

    /// [optional] some plots
    ///
    if(0)
    {
        auto* chart = new TChart();

        /// draw particle interraction force magnitude
        ///
        if(0)
        {
            chart->Scale(Point(10, 0.1));

            auto R = range(0, 2*20, 1000);

            chart->push_back
            (
                new Graph
                (
                    zipmap([](auto r){ return A0(20, r); }, R),
                    pixel(0, 0, 200)
                )
            );
        }

        /// draw strain graph [dynamic]
        ///
        if(1)
        {
            Graph* G[] =
            {
                new Graph(pixel(200, 0, 0)),        /// log of total kinetic energy [red]
                new Graph(pixel(0, 170, 0)),        /// average strain [green]
                new Graph(pixel(200, 100, 0))       /// max strain [orange]
            };

            for(auto* g: G)
                chart->push_back(g);

            *w += [P, G](sEvent& event)
            {
                static double x0 = 0;
                static auto T = clock();
                const double dt = 0.1*CLOCKS_PER_SEC;      /// update time [in ticks]

                switch(event)
                {
                    case sEventCodes::Keyboard_Down:

                        switch(sKey.Key)
                        {
                            case 46:
                                for(auto* g: G)
                                    g->clear();

                                T = clock();
                                x0 = 1.*T/CLOCKS_PER_SEC;
                                break;
                        }
                        break;

                    default:
                        break;
                }

                ///---------------------------------------------------

                auto t = clock();

                if(t - T >= dt)
                {
                    T = t;

                    /// ~ kinetic energy
                    ///
                    double E = average
                    (
                        map([](const auto& p){ return p.v.Sqr(); }, *P)
                    );

                    double S = 0;       /// ~ strain
                    double Z = 0;       /// ~ max strain

                    unsigned N = 0;

                    for(unsigned k = 0; k < P->size(); k++)
                    {
                        for(unsigned p = k+1; p < P->size(); p++)
                        {
                            double dp = (P->at(k).p - P->at(p).p).Norm();
                            double R0 = P->at(k).R0 + P->at(p).R0;

                            if(dp <= R0)
                            {
                                double s = 1 - dp/R0;

                                S += s;
                                N++;

                                if(s > Z)
                                    Z = s;
                            }
                        }
                    }

                    S = 100*S/N;
                    Z = 100*Z;

                    double x = 1.*T/CLOCKS_PER_SEC - x0;

                    G[0]->push_back(Point(x, std::log10(E)));
                    G[1]->push_back(Point(x, S));
                    G[2]->push_back(Point(x, Z));
                }
            };
        }

        *w += chart;
    }

    w->Show();

    return 0;
}
