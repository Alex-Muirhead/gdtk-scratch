/**
 * solidfvcell.d
 *
 * A solid finite-volume cell, to be held by SolidBlock objects.
 *
 * Author: Rowan G. Kyle D. and Peter J.
 * Version: 2015-22-04
 */

module solidfvcell;

import std.conv;
import std.string;
import std.array;
import std.format;
import nm.complex;
import nm.number;
import geom;
import solidfvinterface;
import solidfvvertex;
import solidprops;
import std.stdio;
import globalconfig;

class SolidFVCell {
public:
    size_t id;
    // Cell properties
    number volume;
    number areaxy;
    Vector3 pos;
    // Cell material properties
    SolidProps sp;
    // Cell state
    number T;
    number[] e;
    number[] dedt;
    number de_prev;
    // Cell source term
    number Q;
    // Connections
    SolidFVInterface[] iface;
    SolidFVVertex[] vtx;
    // Cell-centered gradients
    Vector3[] cloud_pos;
    number*[] cloud_T;
    number[12] wx, wy, wz;
    number dTdx;
    number dTdy;
    number dTdz;
    bool is_ghost = true;

private:
    LocalConfig myConfig;

public:

    this(LocalConfig myConfig)
    {
        this.myConfig = myConfig;
        e.length = myConfig.n_flow_time_levels;
        dedt.length = myConfig.n_flow_time_levels;
    }

    @nogc
    void copy_values_from(SolidFVCell other) {
        volume = other.volume;
        areaxy = other.areaxy;
        pos = other.pos;
        sp = other.sp;
        T = other.T;
        foreach (i; 0..e.length) { e[i] = other.e[i] ; }
        foreach (i; 0..dedt.length) { dedt[i] = other.dedt[i] ; }
        de_prev = other.de_prev;
        Q = other.Q;
        dTdx = other.dTdx;
        dTdy = other.dTdy;
        dTdz = other.dTdz;
        is_ghost = other.is_ghost;
    }

    void scanValuesFromString(string buffer)
    {
        auto items = split(buffer);
        pos.x = to!double(items.front); items.popFront();
        pos.y = to!double(items.front); items.popFront();
        pos.z = to!double(items.front); items.popFront();
        volume = to!double(items.front); items.popFront();
        e[0] = to!double(items.front); items.popFront();
        T = to!double(items.front); items.popFront();
        sp.rho = to!double(items.front); items.popFront();
        sp.Cp = to!double(items.front); items.popFront();
        sp.k = to!double(items.front); items.popFront();
        sp.k11 = to!double(items.front); items.popFront();
        sp.k12 = to!double(items.front); items.popFront();
        sp.k13 = to!double(items.front); items.popFront();
        sp.k21 = to!double(items.front); items.popFront();
        sp.k22 = to!double(items.front); items.popFront();
        sp.k23 = to!double(items.front); items.popFront();
        sp.k31 = to!double(items.front); items.popFront();
        sp.k32 = to!double(items.front); items.popFront();
        sp.k33 = to!double(items.front); items.popFront();
    }

    string writeValuesToString() const
    {
        auto writer = appender!string();
        formattedWrite(writer, "%.18e %.18e %.18e %.18e %.18e %.18e %.18e %.18e %.18e %.18e %.18e %.18e %.18e %.18e %.18e %.18e %.18e %.18e",
                       pos.x.re, pos.y.re, pos.z.re, volume.re, e[0].re, T.re,
                       sp.rho.re, sp.Cp.re, sp.k.re,
                       sp.k11.re, sp.k12.re, sp.k13.re,
                       sp.k21.re, sp.k22.re, sp.k23.re,
                       sp.k31.re, sp.k32.re, sp.k33.re);
        return writer.data;
    }


    void timeDerivatives(int ftl, int dimensions)
    {
        SolidFVInterface IFn = iface[Face.north];
        SolidFVInterface IFe = iface[Face.east];
        SolidFVInterface IFs = iface[Face.south];
        SolidFVInterface IFw = iface[Face.west];
        SolidFVInterface IFt, IFb;
         if (dimensions == 3) {
            IFt = iface[Face.top];
            IFb = iface[Face.bottom];
        }
        // Cell volume (inverted).
        number volInv = 1.0 / volume;
        number integral;
        // Sum up fluxes (of form q.n)
        integral = -IFe.flux * IFe.area - IFn.flux * IFn.area
            + IFw.flux * IFw.area + IFs.flux * IFs.area;
        if (dimensions == 3) {
            integral += (-IFt.flux * IFt.area + IFb.flux * IFb.area);
        }
        dedt[ftl] = volInv * integral + Q;
    }

    void stage1RKL1Update(double dt, int j, int s)
    {
        // RKL1 stage 1 coefficients
        double muj_tilde;
        muj_tilde = (2.0*j-1)/j * 2.0/(s*s+s);

        // stage 1 update
        e[1] = e[0] + muj_tilde*dt*dedt[0];

        return;
    } // end stage1RKL1Update()

    void stage2RKL1Update(double dt, int j, int s)
    {
        // stage j coefficients
        double muj; double vuj; double muj_tilde;
        muj_tilde = (2.0*j-1)/j * 2.0/(s*s+s);
        muj = (2.0*j-1)/j;
        vuj = (1.0-j)/j;

        // stage j update
        e[2] = muj*e[1] + vuj*e[0] + muj_tilde*dt*dedt[1];

        return;
    } // end stage2RKL1Update()

    void stage1RKL2Update(double dt, int j, int s)
    {
        // stage 1 coefficients
        double muj_tilde;
        muj_tilde = 4.0/(3.0*(s*s+s-2.0));

        // stage 1 update
        e[1] = e[0] + muj_tilde*dt*dedt[0];

        // make a copy of the initial conserved quantities for the stage j updates
        e[3] = e[0];

        return;
    } // end stage1RKL2Update()

    void stage2RKL2Update(double dt, int j, int s)
    {
        //  stage j coefficients
        double ajm1; double bj; double bjm1, bjm2; double muj; double vuj; double muj_tilde; double gam_tilde;
        if (j == 2) {
            bj = 1.0/3.0;
            bjm1 = 1.0/3.0;
            bjm2 = 1.0/3.0;
        } else if (j == 3) {
            bj = (j*j+j-2.0)/(2.0*j*(j+1.0));
            bjm1 = 1.0/3.0;
            bjm2 = 1.0/3.0;
        } else if (j == 4) {
            bj = (j*j+j-2.0)/(2.0*j*(j+1.0));
            bjm1 = ((j-1.0)*(j-1.0)+(j-1.0)-2.0)/(2.0*(j-1.0)*((j-1.0)+1.0));
            bjm2 = 1.0/3.0;
        } else {
            bj = (j*j+j-2.0)/(2.0*j*(j+1.0));
            bjm1 = ((j-1.0)*(j-1.0)+(j-1.0)-2.0)/(2.0*(j-1.0)*((j-1.0)+1.0));
            bjm2 = ((j-2.0)*(j-2.0)+(j-2.0)-2.0)/(2.0*(j-2.0)*((j-2.0)+1.0));
        }
        ajm1 = 1.0-bjm1;
        muj_tilde = (4.0*(2.0*j-1.0))/(j*(s*s+s-2.0)) * (bj/bjm1);
        gam_tilde = -ajm1*muj_tilde;
        muj = (2*j-1.0)/(j) * (bj/bjm1);
        vuj = -(j-1.0)/(j) * (bj/bjm2);

        // stage j update
        e[2] = muj*e[1] + vuj*e[0] + (1.0-muj-vuj)*e[3] + muj_tilde*dt*dedt[1] + gam_tilde*dt*dedt[0];

        return;
    } // end stage2RKL2Update()

    void eulerUpdate(double dt)
    {
        e[1] = e[0] + dt*dedt[0];
    }

    void stage1Update(double dt)
    {
        double gamma1 = 1.0; // Assume Euler
        final switch (myConfig.gasdynamic_update_scheme) {
        case GasdynamicUpdate.euler:
        case GasdynamicUpdate.backward_euler:
        case GasdynamicUpdate.implicit_rk1:
        case GasdynamicUpdate.moving_grid_1_stage:
        case GasdynamicUpdate.moving_grid_2_stage:
        case GasdynamicUpdate.pc: gamma1 = 1.0; break;
        case GasdynamicUpdate.midpoint: gamma1 = 0.5; break;
        case GasdynamicUpdate.classic_rk3: gamma1 = 0.5; break;
        case GasdynamicUpdate.tvd_rk3: gamma1 = 1.0; break;
        case GasdynamicUpdate.denman_rk3: gamma1 = 8.0/15.0; break;
        case GasdynamicUpdate.rkl1:
        case GasdynamicUpdate.rkl2: assert(false, "invalid option");
        case GasdynamicUpdate.classic_rk4: gamma1 = 0.5; break;
        }
        e[1] = e[0] + dt*gamma1*dedt[0];
    }

    void stage2Update(double dt)
    {
        // Assuming predictor-corrector
        double gamma1 = 0.5;
        double gamma2 = 0.5;
        final switch (myConfig.gasdynamic_update_scheme) {
        case GasdynamicUpdate.euler:
        case GasdynamicUpdate.backward_euler:
        case GasdynamicUpdate.implicit_rk1:
        case GasdynamicUpdate.moving_grid_1_stage: assert(false, "invalid for 1-stage update.");
        case GasdynamicUpdate.moving_grid_2_stage:
        case GasdynamicUpdate.pc: gamma1 = 0.5, gamma2 = 0.5; break;
        case GasdynamicUpdate.midpoint: gamma1 = 0.0; gamma2 = 1.0; break;
        case GasdynamicUpdate.classic_rk3: gamma1 = -1.0; gamma2 = 2.0; break;
        case GasdynamicUpdate.tvd_rk3: gamma1 = 0.25; gamma2 = 0.25; break;
        case GasdynamicUpdate.denman_rk3: gamma1 = -17.0/60.0; gamma2 = 5.0/12.0; break;
        case GasdynamicUpdate.rkl1:
        case GasdynamicUpdate.rkl2: assert(false, "invalid option");
        case GasdynamicUpdate.classic_rk4: gamma1 = 0.0; gamma2 = 0.5; break;
        }
        e[2] = e[0] + dt*(gamma1*dedt[0] + gamma2*dedt[1]);
    }

    void stage3Update(double dt)
    {
        // Assuming TVD_RK3 scheme as done in flow update
        double gamma1 = 1.0/6.0; // presume TVD_RK3 scheme.
        double gamma2 = 1.0/6.0;
        double gamma3 = 4.0/6.0;
        final switch (myConfig.gasdynamic_update_scheme) {
        case GasdynamicUpdate.euler:
        case GasdynamicUpdate.implicit_rk1:
        case GasdynamicUpdate.backward_euler:
        case GasdynamicUpdate.moving_grid_1_stage:
        case GasdynamicUpdate.moving_grid_2_stage:
        case GasdynamicUpdate.pc:
        case GasdynamicUpdate.midpoint:
            assert(false, "invalid for 2-stage update.");
        case GasdynamicUpdate.classic_rk3: gamma1 = 1.0/6.0; gamma2 = 4.0/6.0; gamma3 = 1.0/6.0; break;
        case GasdynamicUpdate.tvd_rk3: gamma1 = 1.0/6.0; gamma2 = 1.0/6.0; gamma3 = 4.0/6.0; break;
            // FIX-ME: Check that we have Andrew Denman's scheme ported correctly.
        case GasdynamicUpdate.denman_rk3: gamma1 = 0.0; gamma2 = -5.0/12.0; gamma3 = 3.0/4.0; break;
        case GasdynamicUpdate.rkl1:
        case GasdynamicUpdate.rkl2: assert(false, "invalid option");
        case GasdynamicUpdate.classic_rk4: gamma1 = 0.0; gamma2 = 0.0; gamma3 = 1.0; break;
        }
        e[3] = e[0] + dt*(gamma1*dedt[0] + gamma2*dedt[1] + gamma3*dedt[2]);
    }

    void stage4Update(double dt)
    {
        double gamma1 = 1.0/6.0; // presume classic_rk4 scheme.
        double gamma2 = 2.0/6.0;
        double gamma3 = 2.0/6.0;
        double gamma4 = 1.0/6.0;
        final switch (myConfig.gasdynamic_update_scheme) {
        case GasdynamicUpdate.euler:
        case GasdynamicUpdate.implicit_rk1:
        case GasdynamicUpdate.backward_euler:
        case GasdynamicUpdate.moving_grid_1_stage:
        case GasdynamicUpdate.moving_grid_2_stage:
        case GasdynamicUpdate.pc:
        case GasdynamicUpdate.midpoint:
        case GasdynamicUpdate.classic_rk3:
        case GasdynamicUpdate.tvd_rk3:
        case GasdynamicUpdate.denman_rk3:
        case GasdynamicUpdate.rkl1:
        case GasdynamicUpdate.rkl2:
        case GasdynamicUpdate.classic_rk4: gamma1 = 1.0/6.0; gamma2 = 2.0/6.0; gamma3 = 2.0/6.0; gamma4 = 1.0/6.0; break;
        }
        e[3] = e[0] + dt*(gamma1*dedt[0] + gamma2*dedt[1] + gamma3*dedt[2] + gamma4*dedt[3]);
    }
    version(complex_numbers) {
        @nogc void clear_imaginary_components()
        // When performing the complex-step Frechet derivative in the Newton-Krylov accelerator,
        // the flowstate values accumulate imaginary components, so we have to start with a clean slate, so to speak.
        {
            e[0].im = 0.0;
            e[1].im = 0.0;
            dedt[0].im = 0.0;
            dedt[1].im = 0.0;
            T.im = 0.0;
            foreach (face; iface) {
                face.T.im = 0.0;
                face.e.im = 0.0;
                face.flux.im = 0.0;
            }
        } // end clear_imaginary_components()
    } // end version(complex)
}

string[] varListForSolidCell()
{
    string[] list;
    list ~= ["pos.x", "pos.y", "pos.z", "volume"];
    list ~= ["e", "T"];
    list ~= ["rho", "Cp", "k"];
    list ~= ["k11", "k12", "k13"];
    list ~= ["k21", "k22", "k23"];
    list ~= ["k31", "k32", "k33"];
    return list;
}
