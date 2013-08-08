/*****************************************************************************
*
* ALPS/betheansatz: Bethe Ansatz calculation using ALPS library
*
* Copyright (C) 2013 by Ryo IGARASHI <rigarash@issp.u-tokyo.ac.jp>
*
* This software is published under the ALPS Application License; you
* can use, redistribute it and/or modify it under the terms of the
* license, either version 1 or (at your option) any later version.
* 
* You should have received a copy of the ALPS Application License
* along with this software; see the file LICENSE. If not, the license
* is also available from http://alps.comp-phys.org/.
*
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
* IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
* FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT 
* SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE 
* FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE, 
* ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
* DEALINGS IN THE SOFTWARE.
*
*****************************************************************************/

#ifndef ALPS_BETHEANSATZ_H_
#define ALPS_BETHEANSATZ_H_

#include <cmath>
#include <cstdint>
#include <iostream>

#include <alps/parameter/parameters.h>
#include <alps/parapack/worker.h>

namespace alps {
namespace betheansatz {

class betheansatz_worker
    : public alps::parapack::abstract_worker
{
public:
    betheansatz_worker(alps::Parameters const& ps)
        : alps::parapack::abstract_worker(),
          params_(ps)
    {}
    void init_observables(alps::Parameters const& /* params */,
                          alps::ObservableSet& /* obs */)
    {}
    bool is_thermalized() const {
        return true;
    }
    double progress() const{
        return 1.0;
    }
    void run(alps::ObservableSet& obs) {
        if (progress() >= 1.0) {
            return;
        }
        double const t = 1.0;
        double const V = params_["V"] / 2.0;
        std::size_t const L = params_["L"];
        std::size_t const N = params_["N"];

        double phi = 0.0;

        double E0  = energy(N,   V, phi);
        double Ep2 = energy(N+2, V, phi);
        double Em2 = energy(N-2, V, phi);

        double uK = (Ep2 + Em2 - 2.0 * E0) * L / (4.0 * std::M_PI);
        std::cout << "u/K = " << uK << std::endl;

        double mu = (Ep2 - Em2) / 4.0;
        std::cout << "mu  = " << mu << std::endl;
    }

private:
    double energy(std::size_t const n, double const V, double const phi) const {
        std::size_t const L = params_["L"];
        std::vector<double> psi(n, 0.0);
        std::vector<double> psi0(n, 0.0);

        double const PI = M_PI;
        double energy  = 0.0;
        double energy0 = 0.0;
        do {
            energy0 = energy;
            energy = 0.0;
            for (int i = 0; i < n; ++i) {
                psi0[i] = 2.0 * PI * n + phi;
                for (int j = 0; j < n; ++j) {
                    double temp = V * std::sin((psi(i) - psi(j)) / 2.0)
                        / (std::cos((psi(i) + psi(j)) / 2.0)
                           + V * std::cos((psi(i) - psi(j)) / 2.0));
                    psi0[i] += 2.0 * std::atan(temp);
                }
                psi0[i] /= L;
            }
            for (int i = 0; i < n; ++i) {
                psi[i] = psi0[i];
                energy -= 2.0 * (V + std::cos(psi0[i]));
            }
            energy += V * L / 2.0;
        } while (std::abs(energy - energy0) < 10e-7);
        return energy;
    }

    alps::Parameters params_;
};

} // namespace betheansatz
} // namespace alps

#endif // ALPS_BETHEANSATZ_H_
