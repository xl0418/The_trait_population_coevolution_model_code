#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <stdexcept>
#include <cmath>
#include <vector>
#include <array>
#include <tuple>
#include <algorithm>
#include <mutex>
#include "rndutils.hpp"


namespace py = pybind11;
using namespace py::literals;


namespace impl
{

  template <typename T>
  using py_array = py::array_t<T, py::array::c_style | py::array::forcecast>;


  template <typename T, typename I>
  py_array<T> make_py_array(I i)
  {
    return py::array_t<T>(static_cast<ssize_t>(i));
  }


  template <typename T, typename I>
  py_array<T> make_py_array(I i, I j)
  {
    return py::array_t<T>(typename py_array<T>::ShapeContainer{ static_cast<ssize_t>(i), static_cast<ssize_t>(j) });
  }


  struct sim_res
  {
    ssize_t valid_time;
    std::vector<double> N;
    std::vector<double> Z;
    std::vector<double> V;
  };


  py::object to_py(sim_res&& C)
  {
    const auto n = static_cast<ssize_t>(C.N.size());
    // NaN-ify
    for (auto i = 0; i < n; ++i)
    {
      if (C.N[i] == 0.0) C.N[i] = C.Z[i] = C.V[i] = Py_NAN;
    }
    return py::dict{
      "sim_time"_a = C.valid_time,
      "N"_a = py_array<double>(n, C.N.data()),
      "Z"_a = py_array<double>(n, C.Z.data()),
      "V"_a = py_array<double>(n, C.V.data())
    };
  }


  py::object to_py(std::vector<sim_res>& R)
  {
    auto py_N = make_py_array<double>(R.size(), R[0].N.size());
    auto py_Z = make_py_array<double>(R.size(), R[0].Z.size());
    auto py_V = make_py_array<double>(R.size(), R[0].V.size());
    auto py_t = make_py_array<ssize_t>(R.size());
    auto pN = py_N.mutable_unchecked();
    auto pZ = py_Z.mutable_unchecked();
    auto pV = py_V.mutable_unchecked();
    auto pt = py_t.mutable_unchecked();
    const auto n = static_cast<ssize_t>(R[0].N.size());
    for (size_t j = 0; j < R.size(); ++j)
    {
      pt(j) = R[j].valid_time;
      if (!R[j].N.empty())
      {
        for (auto i = 0; i < n; ++i)
        {
          if (R[j].N[i] == 0.0) R[j].N[i] = R[j].Z[i] = R[j].V[i] = Py_NAN;
          pN(j, i) = R[j].N[i];
          pZ(j, i) = R[j].Z[i];
          pV(j, i) = R[j].V[i];
        }
      }
      else
      {
        for (auto i = 0; i < n; ++i)
        {
          pN(j, i) = pZ(j, i) = pV(j, i) = Py_NAN;
        }
      }
    }
    return py::dict{
      "sim_time"_a = std::move(py_t),
      "N"_a = std::move(py_N),
      "Z"_a = std::move(py_Z),
      "V"_a = std::move(py_V)
    };
  }


  enum sim_param : int
  {
    gamma = 0,
    a = 1,
    K = 2,
    nu = 3,
    r = 4,
    theta = 5,
    Vmax = 6,
    init_z = 7,
    init_n = 8,
    split_stddev = 9,
    max_param
  };


  template <typename EVENTS>
  sim_res sim_single(const ssize_t T, const ssize_t S, const EVENTS& events, const double* param)
  {
    const auto gamma = param[sim_param::gamma];
    const auto a = param[sim_param::a];
    const auto K = param[sim_param::K];
    const auto nu = param[sim_param::nu];
    const auto r = param[sim_param::r];
    const auto theta = param[sim_param::theta];
    const auto Vmax = param[sim_param::Vmax];

    std::vector<double> Z(S, param[sim_param::init_z]), nZ(S);
    std::vector<double> N(S, 0.0), nN(S);
    std::vector<double> V(S, 1.0 / S), nV(S);
    auto reng = rndutils::make_random_engine_low_entropy<>();
    {
      auto init_n_dist = std::normal_distribution<double>(param[sim_param::init_n]);
      N[0] = init_n_dist(reng);
      N[1] = init_n_dist(reng);
    }
    auto split_dist = std::normal_distribution<double>(0.5, param[sim_param::split_stddev]);
    ssize_t sp = 2;     // existing species
    ssize_t node = 0;
    auto next_event = events.data(node, 0);
    ssize_t sim_time = 0;
    for (; sim_time < T; ++sim_time)
    {
      for (ssize_t i = 0; i < sp; ++i)
      {
        if (N[i] != 0.0)
        {
          double beta = 0.0;
          double sigma = 0.0;
          double sigmasqr = 0.0;
          const double zi = Z[i];
          for (ssize_t j = 0; j < sp; ++j)
          {
            const double zd = zi - Z[j];
            const double t1 = std::exp(-a * zd * zd) * N[j];
            const double t2 = 2.0 * a * zd;
            beta += t1;
            sigma += t2 * t1;
            sigmasqr += t2 * t2 * t1;

          }
          const auto var_trait = V[i] / (2.0 * N[i]);
          const auto dtz = theta - Z[i];
          nZ[i] = Z[i] + V[i] * (2.0 * gamma * dtz + 1.0 / K * sigma) + std::normal_distribution<double>(0.0, var_trait)(reng);
          const auto possion_lambda = N[i] * r * std::exp(-gamma * dtz * dtz + (1.0 - beta / K));
          nN[i] = static_cast<double>(std::poisson_distribution<ssize_t>(possion_lambda)(reng));
          nV[i] = V[i] / 2.0 + 2.0 * N[i] * nu * Vmax / (1.0 + 4.0 * N[i] * nu)
            + V[i] * V[i] * (
              -2.0 * gamma + 4.0 * gamma * gamma * dtz * dtz +
              (1.0 / K) * (sigma - sigmasqr) + 4.0 * gamma / K *
              dtz * sigma + sigma * sigma / (K * K)
              );
          if ((nV[i] < 0.0) || (nN[i] <= 0.0)) 
            return { sim_time, std::move(N), std::move(Z), std::move(V) };
        }
      }
      if ((sim_time + 1) == next_event[0])
      {
        const auto parent = next_event[1];
        const auto daughter = next_event[2];
        if (daughter == -1)
        { // extinction
          const auto ext = next_event[1];
          nN[ext] = N[ext] = 0.0;
        }
        else
        { // speciation
          nZ[daughter] = nZ[parent];
          auto split_ratio = 0.0;
          while ((split_ratio <= 0.0) || (split_ratio >= 1.0))
          {
            split_ratio = split_dist(reng);
          }
          const auto tmp = nN[parent];
          nN[parent] = split_ratio * tmp;
          nN[daughter] = (1.0 - split_ratio) * tmp;
          nV[parent] *= 0.5;
          nV[daughter] = nV[parent];
          ++sp;
        }
        ++node;
        next_event = events.data(node, 0);
      }
      N.swap(nN);
      Z.swap(nZ);
      V.swap(nV);
    }
    return { sim_time, std::move(N), std::move(Z), std::move(V) };
  }


  template <typename EVENTS, typename PARAM>
  py::object sim_multiple(ssize_t T, ssize_t S, const EVENTS& events, const PARAM params)
  {
    const auto sets = params.shape(0);
    auto simres = std::vector<sim_res>(sets);
#pragma omp parallel for
    for (auto s = 0; s < sets; ++s)
    {
      simres[s] = sim_single(T, S, events, params.template unchecked<2>().data(s, 0));
    }
    return to_py(simres);
  }


  py::object sim_dispatch(const py::object& treedata, py_array<double>& params)
  {
    const auto events = py::cast<py_array<ssize_t>>(treedata.attr("events"));
    const auto T = py::cast<ssize_t>(treedata.attr("evo_time"));
    const auto S = py::cast<ssize_t>(treedata.attr("total_species"));
    if ((events.ndim() == 2) && (events.shape(1) == 3))
    {
      if ((params.ndim() == 1) && (params.shape(0) == sim_param::max_param))
      {
        return impl::to_py(sim_single(T, S, events.unchecked<2>(), params.unchecked<1>().data(0)));
      }
      if ((params.ndim() == 2) && (params.shape(1) == sim_param::max_param))
      {
        return sim_multiple(T, S, events.unchecked<2>(), params);
      }
    }
    throw std::runtime_error("parameter don't match");
    return py::none{};
  }


  py_array<ssize_t> discrete_distribution(const py_array<double>& val, ssize_t num)
  {
    using ddist = rndutils::mutable_discrete_distribution<ssize_t, rndutils::all_zero_policy_uni>;
    auto dist = ddist(val.data(), val.data() + val.shape(0));
    auto res = make_py_array<ssize_t>(num);
    auto p = res.mutable_data();
    auto reng = rndutils::make_random_engine_low_entropy<>();
    for (ssize_t i = 0; i < num; ++i, ++p)
    {
      *p = dist(reng);
    }
    return res;
  }


  py_array<double> cauchy_distribution(double a, double b, ssize_t num)
  {
    auto dist = std::cauchy_distribution<>(a, b);
    auto res = make_py_array<double>(num);
    auto p = res.mutable_data();
    auto reng = rndutils::make_random_engine_low_entropy<>();
    for (ssize_t i = 0; i < num; ++i, ++p)
    {
      *p = dist(reng);
    }
    return res;
  }

}


PYBIND11_MODULE(dvtraitsim_cpp, m) {
  m.doc() = R"pbdoc(
        dvtraitsim_cpp plugin
        ---------------------

        .. currentmodule:: dvtraitsim_cpp

        .. autosummary::
           :toctree: _generate

           DVSim
           discrete_distribution
           cauchy_distribution
    )pbdoc";

  m.def("DVSim", &impl::sim_dispatch, "DVSim");
  m.def("discrete_distribution", &impl::discrete_distribution, "discrete_distribution");
  m.def("cauchy_distribution", &impl::cauchy_distribution, "a"_a = 0.0, "b"_a = 1.0, "num"_a, "cauchy_distribution");

  m.attr("all") = 0;
  m.attr("nodes") = 1;
  m.attr("tips") = 2;

#ifdef VERSION_INFO
  m.attr("__version__") = VERSION_INFO;
#else
  m.attr("__version__") = "0.0.1";
#endif
  m.attr("__author__") = "Hanno Hildenbrandt";
}

