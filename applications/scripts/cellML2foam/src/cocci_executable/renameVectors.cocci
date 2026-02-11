// ==========================================================
// Separate STATES and RATES
// ==========================================================

// ---- STATES: NV_Ith_S(y, i) ----
@states@
expression i;
@@
- NV_Ith_S(y, i)
+ STATES[i]

// ---- RATES: NV_Ith_S(ydot, i) ----
@rates@
expression i;
@@
- NV_Ith_S(ydot, i)
+ RATES[i]


// ==========================================================
// REPLACE CONSTANTS
// ==========================================================

@constants@
identifier C =~ "^AC_.*";
@@
- C
+ CONSTANTS[C]


@algebraic@
identifier V =~ "^AV_.*";
@@
- V
+ ALGEBRAIC[V]

