NEURON {
POINT_PROCESS Gap
RANGE vgap
RANGE g, i
NONSPECIFIC_CURRENT i
RANGE synapseID
}
PARAMETER { g = 1 (nanosiemens) }
ASSIGNED {
v (millivolt)
vgap (millivolt)
i (nanoamp)
synapseID
}
BREAKPOINT { i = (v - vgap)*(g*1e-3) }
