# Changelog

## 1.0

- First release: 142 convex QP instances (`rt-qp-0001` – `rt-qp-0142`) from clinical VMAT plans
  (lung, oligometastases, prostate), each with a Warm-IP reference solution. See [docs/qp.md](docs/qp.md).
- `radoptlib` Python package: load, download, select, info, list_instances, evaluate_solution,
  load_solution, to_mps.
- `radoptlib2mps`: converts a folder of HDF5 instances to MPS files; installed with the package and also
  built as a stand-alone program for users without Python. See [mps/README.md](mps/README.md).
