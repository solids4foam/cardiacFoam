# Study inputs and post-processing

`studies/modelMatrix/sweep_scalar_batched.json` is the declarative 2-D
model/backend matrix. It records the myocardial model/tissue rows and the
shared time-integration settings; it contains no scheduler or site-specific
commands.

`postprocess_backend_parity.py` compares two completed outputs. The first
argument is the reference and the second is the candidate. Both must use the
same mesh, stimulus, output times, and tissue timestep.
