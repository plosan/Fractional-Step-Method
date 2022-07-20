
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
