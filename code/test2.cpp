
void computeVelocityU(double* u, const NCMesh m, const double* u_pred, const double* p, const double tstep, const Properties props) {

    int nx = m.getNX();
    int ny = m.getNY();

    for(int i = 1; i < nx; i++)
        for(int j = 1; j < ny+1; j++)
            u[j*(nx+1)+i] = u_pred[j*(nx+1)+i] - tstep / props.rho * (p[j*(nx+2)+i+1] - p[j*(nx+2)+i]) / m.atDistX(i);

}

void computeVelocityV(double* v, const NCMesh m, const double* v_pred, const double* p, const double tstep, const Properties props) {

    int nx = m.getNX();
    int ny = m.getNY();

    for(int i = 1; i < nx+1; i++)
        for(int j = 1; j < ny; j++)
            v[j*(nx+2)+i] = v_pred[j*(nx+2)+i] - tstep / props.rho * (p[(j+1)*(nx+2)+i] - p[j*(nx+2)+i]) / m.atDistY(j);

}


// void computeRuComplex(double* Ru, const NCMesh m, const double* u, const double* v, const Properties props) {
//
//     int nx = m.getNX();
//     int ny = m.getNY();
//
//     // To save the mass flows of the south face
//     // The mass flow through the north face of node (i,j) equals the mass flow through the south face of node (i,j+1)
//     // When the first row is computed (j = 1), my is filled with the north face mass flows. In subsequent rows (j > 1), the mass flow through the
//     // south face of node (i,j) is taken from my[i] and then my[i] is replaced by the mass flow through the north face of node (i,j)
//     double* my = (double*) calloc(nx+1, sizeof(double));
//     if(!my) {
//         printf("Error: could not allocate enough memory for variable my\n");
//         return;
//     }
//
//     // First row (j = 1)
//     int j = 1;
//
//     // Properties that can be updated instead of computing them twice
//     double uW = u[j*(nx+1)];        // West node velocity
//     double uP = u[j*(nx+1)+1];      // Current node velocity
//     double uw = schemeCDS(uW, uP);  // Velocity west face
//     double mw = 0.5 * props.rho * (uW + uP) * m.atSurfX(j); // Mass flow west face
//
//     for(int i = 1; i < nx; i++) {
//         // Nodal velocities
//         double uE = u[j*(nx+1)+i+1];    // East node velocity
//         double uS = u[(j-1)*(nx+1)+i];  // South node velocity
//         double uN = u[(j+1)*(nx+1)+i];  // North node velocity
//
//         // Mass flows
//         double me = 0.5 * props.rho * (uP + uE) * m.atSurfX(j);                                                         // Mass flow east face
//         double mn = props.rho * (v[j*(nx+2)+i] * m.atSemiSurfY(i,1) + v[j*(nx+2)+i+1] * m.atSemiSurfY(i+1,0));          // Mass flow north face
//         double ms = props.rho * (v[(j-1)*(nx+2)+i] * m.atSemiSurfY(i,1) + v[(j-1)*(nx+2)+i+1] * m.atSemiSurfY(i+1,0));  // Mass flow south face
//         my[i] = mn;
//
//         // Face velocities
//         double ue = schemeCDS(uP, uE);
//         double us = uS;                 // Equal to uS because it is the first row
//         double un = schemeCDS(uP, uN);  // This may give problems
//
//         // Operator R(u)
//         double integral1 = -(me * ue - mw * uw + mn * un - ms * us);
//         double integral2 = m.atSurfX(j) * (uE - uP) / m.atDistFaceX(i) - m.atSurfX(j) * (uP - uW) / m.atDistFaceX(i-1);
//         integral2 += m.atSurfY_StaggX(i) * (uN - uP) / m.atDistY(j) - m.atSurfY_StaggX(i) * (uP - uS) / m.atDistY(j-1);
//         integral2 *= props.mu;
//         Ru[j*(nx+1)+i] = (integral1 + integral2) / m.atVolStaggX(i,j);
//
//         // Next node
//         uW = uP;
//         uP = uE;
//         uw = ue;
//         mw = me;
//     }
//
//     // In between rows
//     for(j = 2; j < ny; j++) {
//         // Properties that can be updated instead of computing them twice
//         uW = u[j*(nx+1)];                                   // West node velocity
//         uP = u[j*(nx+1)+1];                                 // Current node velocity
//         uw = schemeCDS(uW, uP);                             // Velocity west face
//         mw = 0.5 * props.rho * (uW + uP) * m.atSurfX(j);    // Mass flow west face
//         // Go through all nodes in the j-th row
//         for(int i = 1; i < nx; i++) {
//             // Nodal velocities
//             double uE = u[j*(nx+1)+i+1];
//             double uS = u[(j-1)*(nx+1)+i];
//             double uN = u[(j+1)*(nx+1)+i];
//
//             // Mass flows
//             double me = 0.5 * props.rho * (uP + uE) * m.atSurfX(j);                                                 // East face mass flow
//             double ms = my[i];                                                                                      // South face mass flow (previously computed as mass flow through north face)
//             double mn = props.rho * (v[j*(nx+2)+i] * m.atSemiSurfY(i,1) + v[j*(nx+2)+i+1] * m.atSemiSurfY(i+1,0));  // North face mass flow
//             my[i] = mn;                                                                                             // Save north face mass flow (to later on be used as mass flow through south face)
//
//             // Face velocities
//             double ue = schemeCDS(uP, uE);
//             double us = schemeCDS(uP, uS);  // This may give problems. It can be optimised by computing it in the previous row
//             double un = schemeCDS(uP, uN);  // This may give problems
//
//             // Operator R(u)
//             double integral1 = -(me * ue - mw * uw + mn * un - ms * us);
//             double integral2 = m.atSurfX(j) * (uE - uP) / m.atDistFaceX(i) - m.atSurfX(j) * (uP - uW) / m.atDistFaceX(i-1);
//             integral2 += m.atSurfY_StaggX(i) * (uN - uP) / m.atDistY(j) - m.atSurfY_StaggX(i) * (uP - uS) / m.atDistY(j-1);
//             integral2 *= props.mu;
//             Ru[j*(nx+1)+i] = (integral1 + integral2) / m.atVolStaggX(i,j);
//
//             // Next node
//             uW = uP;
//             uP = uE;
//             uw = ue;
//             mw = me;
//         }
//     }
//
//     // Last row (j = ny)
//     j = ny;
//
//     uW = u[j*(nx+1)];
//     uP = u[j*(nx+1)+1];
//     uw = schemeCDS(uW, uP);
//     mw = 0.5 * props.rho * (uW + uP) * m.atSurfX(j);
//
//     for(int i = 1; i < nx; i++) {
//         // Nodal velocities
//         double uE = u[j*(nx+1)+i+1];   // East node velocity
//         double uS = u[(j-1)*(nx+1)+i]; // South node velocity
//         double uN = u[(j+1)*(nx+1)+i]; // North node velocity
//
//         // Mass flows
//         double me = 0.5 * props.rho * (uP + uE) * m.atSurfX(j);
//         double ms = my[i];
//         double mn = props.rho * (v[(j+1)*(nx+2)+i] * m.atSemiSurfY(i,1) + v[(j+1)*(nx+2)+i+1] * m.atSemiSurfY(i+1,0));
//
//         // Face velocities
//         double ue = schemeCDS(uP, uE);
//         double us = schemeCDS(uP, uS);  // This may give problems
//         double un = uN;                 // Equal to uN because it is the last row
//
//         // Operator Ru
//         double integral1 = -(me * ue - mw * uw + mn * un - ms * us);
//         double integral2 = m.atSurfX(j) * (uE - uP) / m.atDistFaceX(i) - m.atSurfX(j) * (uP - uW) / m.atDistFaceX(i-1);
//         integral2 += m.atSurfY_StaggX(i) * (uN - uP) / m.atDistY(j) - m.atSurfY_StaggX(i) * (uP - uS) / m.atDistY(j-1);
//         integral2 *= props.mu;
//         Ru[j*(nx+1)+i] = (integral1 + integral2) / m.atVolStaggX(i,j);
//
//         // Next node
//         uW = uP;
//         uP = uE;
//         uw = ue;
//         mw = me;
//     }
// }
//
// void computeRvComplex(double* Rv, const NCMesh m, const double* u, const double* v, const Properties props) {
//
//
//     int nx = m.getNX();
//     int ny = m.getNY();
//
//     double* mx = (double*) calloc(ny+1, sizeof(double));
//     if(!mx) {
//         printf("Error: could not allocate enough memory for variable mx\n");
//         return;
//     }
//
//     // Column i = 1
//     int i = 1;
//     int j = 1;
//     double vS = v[(j-1)*(nx+2)+i];                          // South node y-velocity
//     double vP = v[j*(nx+2)+i];                              // Current node y-velocity
//     double vs = schemeCDS(vP, vS);                          // South face y-velocity
//     double ms = 0.5 * props.rho * (vP + vS) * m.atSurfY(i); // South face mass flow
//
//     for(j = 1; j < ny; j++) {
//
//         // Nodal velocities
//         double vW = v[j*(nx+2)+i-1];    // West node y-velocity
//         double vE = v[j*(nx+2)+i+1];    // East node y-velocity
//         double vN = v[(j+1)*(nx+2)+i];  // North node y-velocity
//
//         // Mass flows
//         double mw = props.rho * (u[j*(nx+1)+i-1] * m.atSemiSurfX(j,1) + u[(j+1)*(nx+1)+i-1] * m.atSemiSurfX(j+1,0));    // West face mass flow
//         double me = props.rho * (u[j*(nx+1)+i] * m.atSemiSurfX(j,1) + u[(j+1)*(nx+1)+i] * m.atSemiSurfX(j+1,0));        // East face mass flow
//         double mn = 0.5 * props.rho * (vP + vN) * m.atSurfY(i);                                                         // North face mass flow
//         mx[j] = me;
//
//         // Face velocities
//         double vw = vW;                 // West face y-velocity
//         double ve = schemeCDS(vP, vE);  // East face y-velocity. This may give problems
//         double vn = schemeCDS(vP, vN);  // North face y-velocity
//
//         // Operator R(v)
//         double integral1 = -(me * ve - mw * vw + mn * vn - ms * vs);
//         double integral2 = m.atSurfX_StaggY(j) * (vE - vP) / m.atDistX(i) - m.atSurfX_StaggY(j) * (vP - vW) / m.atDistX(i-1);
//         integral2 += m.atSurfY(i) * (vN - vP) / m.atDistFaceY(j) - m.atSurfY(i) * (vP - vS) / m.atDistFaceY(j-1);
//         integral2 *= props.mu;
//         Rv[j*(nx+2)+i] = (integral1 + integral2) / m.atVolStaggY(i,j);
//
//         // Next node
//         vS = vP;
//         vP = vE;
//         vs = vn;
//         ms = mn;
//
//     }
//
//     // Inner columns (2 <= i < nx)
//     for(i = 2; i < nx; i++) {
//
//         vS = v[(j-1)*(nx+2)+i];                             // South node y-velocity
//         vP = v[j*(nx+2)+i];                                 // Current node y-velocity
//         vs = schemeCDS(vP, vS);                             // South face y-velocity
//         ms = 0.5 * props.rho * (vP + vS) * m.atSurfY(i);    // South face mass flow
//
//         for(int j = 1; j < ny; j++) {
//
//             // Nodal velocities
//             double vW = v[j*(nx+2)+i-1];
//             double vE = v[j*(nx+2)+i+1];
//             double vN = v[(j+1)*(nx+2)+i];
//
//             // Mass flows
//             double mw = mx[j];                                                                                          // West face mass flow
//             double me = props.rho * (u[j*(nx+1)+i] * m.atSemiSurfX(j,1) + u[(j+1)*(nx+1)+i] * m.atSemiSurfX(j+1,0));    // West face mass flow
//             double mn = 0.5 * props.rho * (vP + vN) * m.atSurfY(i);                                                     // North face mass flow
//             mx[j] = me;
//
//             // Face velocities
//             double vw = schemeCDS(vP, vW);  // West face y-velocity. This may give problems
//             double ve = schemeCDS(vP, vE);  // East face y-velocity. This may give problems
//             double vn = schemeCDS(vP, vN);  // North face y-velocity. This may give problems
//
//             // Operator R(v)
//             double integral1 = -(me * ve - mw * vw + mn * vn - ms * vs);
//             double integral2 = m.atSurfX_StaggY(j) * (vE - vP) / m.atDistX(i) - m.atSurfX_StaggY(j) * (vP - vW) / m.atDistX(i-1);
//             integral2 += m.atSurfY(i) * (vN - vP) / m.atDistFaceY(j) - m.atSurfY(i) * (vP - vS) / m.atDistFaceY(j-1);
//             integral2 *= props.mu;
//             Rv[j*(nx+2)+i] = (integral1 + integral2) / m.atVolStaggY(i,j);
//
//             // Next node
//             vS = vP;
//             vP = vN;
//             vs = vn;
//             ms = mn;
//
//         }
//
//
//     }
//
//     // Last column (i = nx)
//     i = nx;
//
//     for(j = 1; j < ny; j++) {
//
//         // Nodal velocities
//         double vW = v[j*(nx+2)+i-1];    // West node y-velocity
//         double vE = v[j*(nx+2)+i+1];    // East node y-velocity
//         double vN = v[(j+1)*(nx+2)+i];  // North node y-velocity
//
//         // Mass flows
//         double mw = mx[j];
//         double me = props.rho * (u[j*(nx+1)+i] * m.atSemiSurfX(j,1) + u[(j+1)*(nx+1)+i] * m.atSemiSurfX(j+1,0));
//         double mn = 0.5 * props.rho * (vP + vN) * m.atSurfY(i);
//
//         // Face velocities
//         double vw = schemeCDS(vP, vW);  // West face y-velocity. This may give problems
//         double ve = schemeCDS(vP, vE);  // East face y-velocity. This may give problems
//         double vn = schemeCDS(vP, vN);  // North face y-velocity
//
//         // Operator R(v)
//         double integral1 = -(me * ve - mw * vw + mn * vn - ms * vs);
//         double integral2 = m.atSurfX_StaggY(j) * (vE - vP) / m.atDistX(i) - m.atSurfX_StaggY(j) * (vP - vW) / m.atDistX(i-1);
//         integral2 += m.atSurfY(i) * (vN - vP) / m.atDistFaceY(j) - m.atSurfY(i) * (vP - vS) / m.atDistFaceY(j-1);
//         integral2 *= props.mu;
//         Rv[j*(nx+2)+i] = (integral1 + integral2) / m.atVolStaggY(i,j);
//
//         // Next node
//         vS = vP;
//         vP = vN;
//         vs = vn;
//         ms = mn;
//     }
// }
