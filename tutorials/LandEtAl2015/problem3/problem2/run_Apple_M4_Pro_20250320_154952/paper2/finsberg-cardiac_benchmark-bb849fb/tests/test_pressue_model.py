import numpy as np

from cardiac_benchmark.pressure_model import default_parameters_benchmark1
from cardiac_benchmark.pressure_model import pressure_function
from cardiac_benchmark.problem import LVProblem


def test_active_model_default_parameters():
    problem_parameters = LVProblem.default_parameters()
    dt = float(problem_parameters["dt"])
    time = np.arange(dt, 1, dt)
    t_eval = time - float(problem_parameters["alpha_f"]) * dt
    pressure = pressure_function(
        t_span=(0, 1),
        t_eval=t_eval,
        parameters=default_parameters_benchmark1(),
    )
    # Different value than in the paper, but it depends on dt, alpha_f and integration method
    assert np.isclose(pressure.max(), 16073.142358747114)

    import matplotlib.pyplot as plt

    plt.plot(time, pressure)
    plt.savefig("pressure.png")
