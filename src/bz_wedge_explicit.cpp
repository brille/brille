/*
  This method is a strongly-modified version of the function KPOINTS_Directions
  from AFLOW 3.1.244 -- http://materials.duke.edu/awrapper.html
  https://arxiv.org/abs/1004.2974 -- which was licences under GPL v3.
*/
bool BrillouinZone::wedge_explicit(void){
  Reciprocal rlat = this->outerlattice;
  Direct dlat = rlat.star();
  Spacegroup spgr = rlat.get_spacegroup_object();
  Pointgroup ptgr = rlat.get_pointgroup_object();
  Holohedry holo = ptgr.get_holohedry();
  Bravais brav = spgr.bravais;

  double alpha = dlat.get_alpha(), beta=dlat.get_beta(), gamma=dlat.get_gamma();
  double a = dlat.get_a(), b=dlat.get_b(), c = dlat.get_c();

  double c_alpha = std::cos(alpha), s_alpha = std::sin(alpha);
  double c_beta = std::cos(beta), s_beta = std::sin(beta);
  double c_gamma = std::cos(gamma), s_gamma = std::sin(gamma);
  double t_alpha = s_alpha/c_alpha, t_beta = s_beta/c_beta, t_gamma = s_gamma/c_gamma;

  double ap, h2, h, ka, kalpha, kap, kh2, kh;
  double zeta, eta, mu, nu, dlt, psi, phi, omega, rho;

  std::vector<std::array<double, 3>> points;
  points.push_back({0.,0.,0.});
  switch (holo){
    case Holohedry::cubic:
      if (brav == Bravais::P){
        points.push_back({0.500, 0.500, 0.000});
        points.push_back({0.500, 0.500, 0.500});
        points.push_back({0.000, 0.500, 0.000});
      }
      if (brav == Bravais::F){
        points.push_back({0.375, 0.375, 0.750});
        points.push_back({0.500, 0.500, 0.500});
        points.push_back({0.625, 0.250, 0.625});
        points.push_back({0.500, 0.250, 0.750});
        points.push_back({0.500, 0.000, 0.500});
      }
      if (brav == Bravais::I){
        points.push_back({0.000, 0.000, 0.000});
        points.push_back({0.500,-0.500, 0.500});
        points.push_back({0.000, 0.000, 0.500});
        points.push_back({0.250, 0.250, 0.250});
      }
      if (brav != Bravais::P && brav != Bravais::F && brav != Bravais::I)
          throw std::logic_error("Unknown cubic centring");
      break;
    case Holohedry::hexagonal:
      points.push_back({0.000, 0.000, 0.000});
      points.push_back({0.000, 0.000, 0.500});
      points.push_back({1./3., 1./3., 0.500});
      points.push_back({1./3., 1./3., 0.000});
      points.push_back({0.500, 0.000, 0.500});
      points.push_back({0.500, 0.000, 0.000});
      break;
    case Holohedry::trigonal:
      if (!spgr.get_choice().compare("R")){ // Rhombohedral
        ap = 2.*a*std::sin(alpha/2.);
        h2 = a*std::cos(alpha/2.);
        h2 *= h2; // poor-man's squaring of a complicated number
        h2 -= ap*ap/12;
        h = std::sqrt(h2);
        if (2*alpha < PI){
          eta = 5./6. - ap*ap/h2/18.;
          nu = (3./2.-eta)/2.;
          points.push_back({0.500, 0.500, 0.000});
          points.push_back({0.500, 0.000, 0.000});
          points.push_back({0.000, 0.000,-0.500});
          points.push_back({0.500, 0.500, 0.500});
          points.push_back({  eta, 0.500, 1-eta});
          points.push_back({0.500, 1-eta, eta-1});
          points.push_back({  eta,    nu,    nu});
          points.push_back({ 1-nu,  1-nu, 1-eta});
          points.push_back({   nu,    nu, eta-1});
          points.push_back({ 1-nu,    nu, 0.000});
          points.push_back({   nu, 0.000,   -nu});
        } else {
          ka = std::sqrt(4./3/ap/ap+1./9/h2);
          kalpha = std::acos((1.0/9/h2-2.0/3/ap/ap)/ka/ka);
          kap = 2.*ka*std::sin(kalpha/2.);
          kh2 = ka*std::cos(kalpha/2.);
          kh2 *= kh2;
          kh2 -= kap*kap/12.;
          kh = std::sqrt(kh2);
          eta = 0.5*ka*ka*h/kh;
          nu = (3./2.-eta)/2.;
          points.push_back({0.500,-0.500, 0.000});
          points.push_back({0.500, 0.000, 0.000});
          points.push_back({0.500,-0.500, 0.500});
          points.push_back({  eta,   eta,   eta});
          points.push_back({1-eta,  -eta,  -eta});
          points.push_back({ 1-nu,   -nu,  1-nu});
          points.push_back({   nu,  nu-1,  nu-1});
        }
      } else { // Hexagonal
        points.push_back({0.000, 0.000, 0.000});
        points.push_back({0.000, 0.000, 0.500});
        points.push_back({1./3., 1./3., 0.500});
        points.push_back({1./3., 1./3., 0.000});
        points.push_back({0.500, 0.000, 0.500});
        points.push_back({0.500, 0.000, 0.000});
      }
      break;
    case Holohedry::tetragonal:
      if (brav == Bravais::P){
        points.push_back({0.500, 0.500, 0.500});
        points.push_back({0.500, 0.500, 0.000});
        points.push_back({0.000, 0.500, 0.500});
        points.push_back({0.000, 0.500, 0.000});
        points.push_back({0.000, 0.000, 0.500});
      }
      if (brav == Bravais::I){
        if (a < c){
          eta = (1+a*a/c/c)/4.;
          zeta = a*a/c/c/2.;
          points.push_back({0.250, 0.250, 0.250});
          points.push_back({0.000, 0.000, 0.500});
          points.push_back({0.500, 0.500,-0.500});
          points.push_back({ -eta,   eta,   eta});
          points.push_back({  eta, 1-eta,  -eta});
          points.push_back({-zeta,  zeta, 0.500});
          points.push_back({0.500, 0.500, -zeta});
        } else {
          eta = (1+c*c/a/a)/4.;
          points.push_back({-0.500,0.500, 0.500});
          points.push_back({0.000, 0.500, 0.000});
          points.push_back({0.250, 0.250, 0.250});
          points.push_back({0.000, 0.000, 0.500});
          points.push_back({  eta,   eta,  -eta});
          points.push_back({ -eta, 1-eta,   eta});
        }
      }
      if (brav != Bravais::I && brav != Bravais::P)
          throw std::logic_error("Unknown tetragonal centring");
    break;
    case Holohedry::orthogonal:
      if (brav == Bravais::P){
        points.push_back({0.500, 0.500, 0.500});
        points.push_back({0.500, 0.500, 0.000});
        points.push_back({0.000, 0.500, 0.500});
        points.push_back({0.500, 0.000, 0.500});
        points.push_back({0.500, 0.000, 0.000});
        points.push_back({0.000, 0.500, 0.000});
        points.push_back({0.000, 0.000, 0.500});
      }
      if (brav == Bravais::F){
        if ( 1/a/a < (1/b/b + 1/c/c) ){
          eta = 0.25*(1+a*a/b/b-a*a/c/c);
          phi = 0.25*(1+c*c/b/b-c*c/a/a);
          dlt = 0.25*(1+b*b/a/a-b*b/c/c);
          points.push_back({  0.500,   0.500, 0.500});
          points.push_back({  0.000,   0.500, 0.500});
          points.push_back({  0.500,   0.000, 0.500});
          points.push_back({  0.500,   0.500, 0.000});
          points.push_back({  0.500, 0.5-eta, 1-eta});
          points.push_back({  0.500, 0.5+eta,   eta});
          points.push_back({0.5-dlt,   0.500, 1-dlt});
          points.push_back({0.5+dlt,   0.500,   dlt});
          points.push_back({  1-phi, 0.5-phi, 0.500});
          points.push_back({    phi, phi+0.5, 0.500});
        } else {
          zeta = 0.25*(1+a*a/b/b-a*a/c/c);
          eta = 0.25*(1+a*a/b/b+a*a/c/c);
          points.push_back({0.500,    0.500,  0.500});
          points.push_back({1.000,    0.500,  0.500});
          points.push_back({0.500,    0.000,  0.500});
          points.push_back({0.500,    0.500,  0.000});
          points.push_back({0.500, 0.5+zeta,   zeta});
          points.push_back({0.500, 0.5-zeta, 1-zeta});
          points.push_back({0.000,      eta,    eta});
          points.push_back({1.000,    1-eta,  1-eta});
        }
      }
      if (brav == Bravais::I){
        zeta= 0.25*(1+a*a/c/c);
        eta = 0.25*(1+b*b/c/c);
        dlt = 0.25*(b*b-a*a)/c/c;
        mu  = 0.25*(a*a+b*b)/c/c;
        points.push_back({  0.000,   0.500,   0.000});
        points.push_back({  0.500,   0.000,   0.000});
        points.push_back({  0.000,   0.000,   0.500});
        points.push_back({  0.250,   0.250,   0.250});
        points.push_back({  0.500,   0.500,-  0.500});
        points.push_back({    -mu,      mu, 0.5-dlt});
        points.push_back({     mu,     -mu, 0.5+dlt});
        points.push_back({0.5-dlt, 0.5+dlt,     -mu});
        points.push_back({  -zeta,    zeta,    zeta});
        points.push_back({   zeta,  1-zeta,   -zeta});
        points.push_back({    eta,    -eta,     eta});
        points.push_back({  1-eta,     eta,    -eta});
      }
      if (brav == Bravais::C){
        if (a<b){
          zeta=(a*a+b*b)/(4.0*b*b);
          points.push_back({ 0.000, 0.500, 0.500});
          points.push_back({ 0.000, 0.500, 0.000});
          points.push_back({-0.500, 0.500, 0.500});
          points.push_back({-0.500, 0.500, 0.000});
          points.push_back({ 0.000, 0.000, 0.500});
          points.push_back({  zeta,  zeta, 0.500});
          points.push_back({ -zeta,1-zeta, 0.500});
          points.push_back({  zeta,  zeta, 0.000});
          points.push_back({ -zeta,1-zeta, 0.000});
        } else {
          zeta=(a*a+b*b)/(4.0*a*a);
          points.push_back({ 0.500, 0.000, 0.500});
          points.push_back({ 0.500, 0.000, 0.000});
          points.push_back({ 0.500,-0.500, 0.500});
          points.push_back({ 0.500,-0.500, 0.000});
          points.push_back({ 0.000, 0.000, 0.500});
          points.push_back({  zeta,  zeta, 0.500});
          points.push_back({1-zeta, -zeta, 0.500});
          points.push_back({  zeta,  zeta, 0.000});
          points.push_back({1-zeta, -zeta, 0.000});
        }
      }
      if (brav == Bravais::A){
        if (b<c){
          zeta=(b*b+c*c)/(4.0*c*c);
          points.push_back({ 0.500, 0.000, 0.500});
          points.push_back({ 0.000, 0.000, 0.500});
          points.push_back({ 0.500,-0.500, 0.500});
          points.push_back({ 0.000,-0.500, 0.500});
          points.push_back({ 0.500, 0.000, 0.000});
          points.push_back({ 0.500,  zeta,  zeta});
          points.push_back({ 0.500, -zeta,1-zeta});
          points.push_back({ 0.000,  zeta,  zeta});
          points.push_back({ 0.000, -zeta,1-zeta});
        } else {
          zeta=(b*b+c*c)/(4.0*b*b);
          points.push_back({ 0.500, 0.500, 0.000});
          points.push_back({ 0.000, 0.500, 0.000});
          points.push_back({ 0.500, 0.500,-0.500});
          points.push_back({ 0.000, 0.500,-0.500});
          points.push_back({ 0.500, 0.000, 0.000});
          points.push_back({ 0.500,  zeta,  zeta});
          points.push_back({ 0.500,1-zeta, -zeta});
          points.push_back({ 0.000,  zeta,  zeta});
          points.push_back({ 0.000,1-zeta, -zeta});
        }
      }
      if (brav == Bravais::B){
        if (c<a){
          zeta=(c*c+a*a)/(4.0*a*a);
          points.push_back({ 0.000, 0.500, 0.500});
          points.push_back({ 0.000, 0.000, 0.500});
          points.push_back({-0.500, 0.500, 0.500});
          points.push_back({-0.500, 0.000, 0.500});
          points.push_back({ 0.000, 0.500, 0.000});
          points.push_back({  zeta, 0.500,  zeta});
          points.push_back({ -zeta, 0.500,1-zeta});
          points.push_back({  zeta, 0.000,  zeta});
          points.push_back({ -zeta, 0.000,1-zeta});
        } else {
          zeta=(c*c+a*a)/(4.0*c*c);
          points.push_back({ 0.500, 0.500, 0.000});
          points.push_back({ 0.500, 0.000, 0.000});
          points.push_back({ 0.500, 0.500,-0.500});
          points.push_back({ 0.500, 0.000,-0.500});
          points.push_back({ 0.000, 0.500, 0.000});
          points.push_back({  zeta, 0.500,  zeta});
          points.push_back({1-zeta, 0.500, -zeta});
          points.push_back({  zeta, 0.000,  zeta});
          points.push_back({1-zeta, 0.000, -zeta});
        }
      }
      if (brav == Bravais::_)
          throw std::logic_error("Unknown orthorhombic centring");
    break;
    case Holohedry::monoclinic:
      if (brav == Bravais::I || brav == Bravais::P){
        if (spgr.get_choice().find("a") != std::string::npos){
          eta = 0.5*(1-b*c_alpha/c)/s_alpha/s_alpha;
          nu = 0.5-eta*c*c_alpha/b;
          points.push_back({ 0.500, 0.500, 0.000});
          points.push_back({ 0.000, 0.500, 0.500});
          points.push_back({ 0.500, 0.000, 0.500});
          points.push_back({ 0.500, 0.000,-0.500});
          points.push_back({ 0.500, 0.500, 0.500});
          points.push_back({ 0.000, 0.500, 0.000});
          points.push_back({ 0.000, 0.000, 0.500});
          points.push_back({ 0.000, 0.000,-0.500});
          points.push_back({ 0.500, 0.000, 0.000});
          points.push_back({ 0.000,   eta, 1.-nu});
          points.push_back({ 0.000, 1-eta,    nu});
          points.push_back({ 0.000,   eta,   -nu});
          points.push_back({ 0.500,   eta, 1.-nu});
          points.push_back({ 0.500, 1-eta,    nu});
          points.push_back({ 0.500,   eta,   -nu});
        }
        if (spgr.get_choice().find("b") != std::string::npos){
          eta = 0.5*(1-c*c_beta/a)/s_beta/s_beta;
          nu = 0.5-eta*a*c_beta/c;
          points.push_back({ 0.500, 0.500, 0.000});
          points.push_back({ 0.500, 0.000, 0.500});
          points.push_back({ 0.000, 0.500, 0.500});
          points.push_back({ 0.000, 0.500,-0.500});
          points.push_back({ 0.500, 0.500, 0.500});
          points.push_back({ 0.500, 0.000, 0.000});
          points.push_back({ 0.000, 0.000, 0.500});
          points.push_back({ 0.000, 0.000,-0.500});
          points.push_back({ 0.000, 0.500, 0.000});
          points.push_back({   eta, 0.000, 1.-nu});
          points.push_back({ 1-eta, 0.000,    nu});
          points.push_back({   eta, 0.000,   -nu});
          points.push_back({   eta, 0.500, 1.-nu});
          points.push_back({ 1-eta, 0.500,    nu});
          points.push_back({   eta, 0.500,   -nu});
        }
        if (spgr.get_choice().find("c") != std::string::npos){
          eta = 0.5*(1-a*c_gamma/b)/s_gamma/s_gamma;
          nu = 0.5-eta*b*c_gamma/a;
          points.push_back({ 0.500, 0.000, 0.500});
          points.push_back({ 0.500, 0.500, 0.000});
          points.push_back({ 0.000, 0.500, 0.500});
          points.push_back({ 0.000,-0.500, 0.500});
          points.push_back({ 0.500, 0.500, 0.500});
          points.push_back({ 0.500, 0.000, 0.000});
          points.push_back({ 0.000, 0.500, 0.000});
          points.push_back({ 0.000,-0.500, 0.000});
          points.push_back({ 0.000, 0.000, 0.500});
          points.push_back({   eta, 1.-nu, 0.000});
          points.push_back({ 1-eta,    nu, 0.000});
          points.push_back({   eta,   -nu, 0.000});
          points.push_back({   eta, 1.-nu, 0.500});
          points.push_back({ 1-eta,    nu, 0.500});
          points.push_back({   eta,   -nu, 0.500});
        }
      }
      if (brav == Bravais::C){
        if (spgr.get_choice().find("a") != std::string::npos){
          if (2*rlat.get_gamma() < PI){
            if (c_alpha*b/c + s_alpha*s_alpha*b*b/a/a > 1){
              zeta = 0.25/s_alpha/s_alpha + 0.25*b*b/a/a - 0.25*b/c/s_alpha/t_alpha;
              eta = 0.5 + 2.0*c/b*zeta*c_alpha;
              mu = 0.25 + 0.25*b*b/a/a - 0.25/t_alpha/t_alpha + 0.25*c*c_alpha/b*(1.0/s_alpha/s_alpha-b*b/a/a);
              phi = 0.75 - 0.25*b*b/a/a - 0.25*b/t_alpha*(1.0/c/s_alpha-1.0/b/t_alpha);
              psi = 0.5 + (2.0*phi-1.0)*c/b*c_alpha;
              nu = 0.25 + 0.25*(b*b/a/a+b/c/s_alpha/t_alpha-3/t_alpha/t_alpha) + 0.5*c*c_alpha/b*(1/s_alpha/s_alpha-b*b/a/a);
              omega = c*c_alpha/b*(2.0*nu-0.5) + c*s_alpha*t_alpha/b*(2.0*nu-0.5*b*b/a/a-0.5);
              rho = 0.75 + 0.25*a*a/b/b/s_alpha*(b/c/t_alpha-1/s_alpha);
              dlt = (2.0*mu-nu)*c*c_alpha/b + 0.5*omega - 0.25;
              points.push_back({  0.500,  0.500,  0.500});
              points.push_back({  0.500,  0.000,  0.500});
              points.push_back({  0.500,  0.000,  0.000});
              points.push_back({  0.000, -0.500,  0.000});
              points.push_back({  0.500, -0.500,  0.000});
              points.push_back({  0.000,  0.000,  0.500});
              points.push_back({     nu,     nu,  omega});
              points.push_back({  1.-nu,  1.-nu,1-omega});
              points.push_back({     nu,  nu-1.,  omega});
              points.push_back({   zeta,   zeta,    eta});
              points.push_back({1.-zeta,  -zeta, 1.-eta});
              points.push_back({  -zeta,  -zeta, 1.-eta});
              points.push_back({    rho, 1.-rho,  0.500});
              points.push_back({ 1.-rho, rho-1.,  0.500});
              points.push_back({    phi, phi-1.,    psi});
              points.push_back({ 1.-phi, 1.-phi, 1.-psi});
              points.push_back({ 1.-phi,   -phi, 1.-psi});
              points.push_back({     mu,     mu,    dlt});
              points.push_back({  1.-mu,    -mu,   -dlt});
              points.push_back({    -mu,    -mu,   -dlt});
              points.push_back({     mu,   mu-1,    dlt});
            } else {
              zeta = 0.25/s_alpha/s_alpha + 0.25*b*b/a/a - 0.25*b/c/s_alpha/t_alpha;
              eta = 0.5 + 2*c/b*zeta*c_alpha;
              mu = 0.25 + 0.25*b*b/a/a;
              dlt = (2*mu-0.5)*c/b*c_alpha;
              phi = 0.75 - 0.25*b*b/a/a - 0.25*b/t_alpha*(1/c/s_alpha-1/b/t_alpha);
              psi = 0.5 + (2*phi-1)*c/b*c_alpha;
              points.push_back({  0.500, -0.500,  0.500});
              points.push_back({  0.500,  0.000,  0.500});
              points.push_back({  0.500,  0.000,  0.000});
              points.push_back({  0.000, -0.500,  0.000});
              points.push_back({  0.500, -0.500,  0.000});
              points.push_back({  0.000,  0.000,  0.500});
              points.push_back({ 1.-phi, 1.-phi, 1.-psi});
              points.push_back({    phi, phi-1.,    psi});
              points.push_back({ 1.-phi,   -phi, 1.-psi});
              points.push_back({   zeta,   zeta,    eta});
              points.push_back({1.-zeta,  -zeta, 1.-eta});
              points.push_back({  -zeta,  -zeta, 1.-eta});
              points.push_back({     mu,     mu,    dlt});
              points.push_back({  1.-mu,    -mu,   -dlt});
              points.push_back({    -mu,    -mu,   -dlt});
              points.push_back({     mu,  mu-1.,    dlt});
            }
          } else {
            zeta = 0.5/s_alpha/s_alpha - 0.25*b/c/t_alpha/s_alpha;
            eta = 0.5 + 2*c*c_alpha*zeta/b;
            psi = 0.75 - 0.25*a*a/b/b/s_alpha/s_alpha;
            phi = psi + 0.25*a*a/b/c/t_alpha/s_alpha;
            mu = psi + 0.25*(2-b/c/c_alpha)/t_alpha/t_alpha;
            dlt = 1.0 - 0.5*b/t_alpha*(1/c/s_alpha-2/b/t_alpha) - mu;
            points.push_back({  0.000,  0.000,  0.000});
            points.push_back({  0.500,  0.500,  0.500});
            points.push_back({  0.500,  0.000,  0.500});
            points.push_back({  0.500,  0.000,  0.000});
            points.push_back({  0.000, -0.500,  0.000});
            points.push_back({  0.500,  0.500,  0.000});
            points.push_back({ -0.500, -0.500,  0.000});
            points.push_back({  0.000,  0.000,  0.500});
            points.push_back({1.-zeta,1.-zeta, 1.-eta});
            points.push_back({   zeta,   zeta,    eta});
            points.push_back({  -zeta,  -zeta, 1.-eta});
            points.push_back({ 1.-dlt,  1.-mu, 1.-eta});
            points.push_back({   -dlt,    -mu, 1.-eta});
            points.push_back({    phi, 1.-phi,  0.500});
            points.push_back({ 1.-phi, phi-1.,  0.500});
            points.push_back({  1.-mu,   -dlt, 1.-eta});
            points.push_back({     mu,    dlt,    eta});
            points.push_back({ 1.-psi, psi-1.,  0.000});
            points.push_back({    psi, 1.-psi,  0.000});
            points.push_back({  psi-1.,  -psi,  0.000});
          }
        }
        if (spgr.get_choice().find("b") != std::string::npos){
          if (2*rlat.get_gamma() < PI){
            if (c_beta*c/a + s_beta*s_beta*c*c/b/b > 1){
              zeta = 0.25/s_beta/s_beta + 0.25*c*c/b/b - 0.25*c/a/s_beta/t_beta;
              eta = 0.5 + 2.0*a/c*zeta*c_beta;
              mu = 0.25 + 0.25*c*c/b/b - 0.25/t_beta/t_beta + 0.25*a*c_beta/c*(1.0/s_beta/s_beta-c*c/b/b);
              phi = 0.75 - 0.25*c*c/b/b - 0.25*c/t_beta*(1.0/a/s_beta-1.0/c/t_beta);
              psi = 0.5 + (2.0*phi-1.0)*a/c*c_beta;
              nu = 0.25 + 0.25*(c*c/b/b+c/a/s_beta/t_beta-3/t_beta/t_beta) + 0.5*a*c_beta/c*(1/s_beta/s_beta-c*c/b/b);
              omega = a*c_beta/c*(2.0*nu-0.5) + a*s_beta*t_beta/c*(2.0*nu-0.5*c*c/b/b-0.5);
              rho = 0.75 + 0.25*b*b/c/c/s_beta*(c/a/t_beta-1/s_beta);
              dlt = (2.0*mu-nu)*a*c_beta/c + 0.5*omega - 0.25;
              points.push_back({ 0.500,  0.500,  0.500});
              points.push_back({ 0.000,  0.500,  0.500});
              points.push_back({ 0.000,  0.000,  0.500});
              points.push_back({-0.500,  0.000,  0.000});
              points.push_back({-0.500,  0.000,  0.500});
              points.push_back({ 0.000,  0.500,  0.000});
              points.push_back({    nu,  omega,     nu});
              points.push_back({ 1.-nu,1-omega,  1.-nu});
              points.push_back({ nu-1.,  omega,     nu});
              points.push_back({  zeta,    eta,   zeta});
              points.push_back({ -zeta, 1.-eta,1.-zeta});
              points.push_back({ -zeta, 1.-eta,  -zeta});
              points.push_back({1.-rho,  0.500,    rho});
              points.push_back({rho-1.,  0.500, 1.-rho});
              points.push_back({phi-1.,    psi,    phi});
              points.push_back({1.-phi, 1.-psi, 1.-phi});
              points.push_back({  -phi, 1.-psi, 1.-phi});
              points.push_back({    mu,    dlt,     mu});
              points.push_back({   -mu,   -dlt,  1.-mu});
              points.push_back({   -mu,   -dlt,    -mu});
              points.push_back({  mu-1,    dlt,     mu});
            } else {
              zeta = 0.25/s_beta/s_beta + 0.25*c*c/b/b - 0.25*c/a/s_beta/t_beta;
              eta = 0.5 + 2*a/c*zeta*c_beta;
              mu = 0.25 + 0.25*c*c/b/b;
              dlt = (2*mu-0.5)*a/c*c_beta;
              phi = 0.75 - 0.25*c*c/b/b - 0.25*c/t_beta*(1/a/s_beta-1/c/t_beta);
              psi = 0.5 + (2*phi-1)*a/c*c_beta;
              points.push_back({-0.500,  0.500,  0.500});
              points.push_back({ 0.000,  0.500,  0.500});
              points.push_back({ 0.000,  0.000,  0.500});
              points.push_back({-0.500,  0.000,  0.000});
              points.push_back({-0.500,  0.000,  0.500});
              points.push_back({ 0.000,  0.500,  0.000});
              points.push_back({1.-phi, 1.-psi, 1.-phi});
              points.push_back({phi-1.,    psi,    phi});
              points.push_back({  -phi, 1.-psi, 1.-phi});
              points.push_back({  zeta,    eta,   zeta});
              points.push_back({ -zeta, 1.-eta,1.-zeta});
              points.push_back({ -zeta, 1.-eta,  -zeta});
              points.push_back({    mu,    dlt,     mu});
              points.push_back({   -mu,   -dlt,  1.-mu});
              points.push_back({   -mu,   -dlt,    -mu});
              points.push_back({ mu-1.,    dlt,     mu});
            }
          } else {
            zeta = 0.5/s_beta/s_beta - 0.25*c/a/t_beta/s_beta;
            eta = 0.5 + 2*a*c_beta*zeta/c;
            psi = 0.75 - 0.25*b*b/c/c/s_beta/s_beta;
            phi = psi + 0.25*b*b/c/a/t_beta/s_beta;
            mu = psi + 0.25*(2-c/a/c_beta)/t_beta/t_beta;
            dlt = 1.0 - 0.5*c/t_beta*(1/a/s_beta-2/c/t_beta) - mu;
            points.push_back({  0.000,  0.000,   0.000});
            points.push_back({  0.500,  0.500,   0.500});
            points.push_back({  0.000,  0.500,   0.500});
            points.push_back({  0.000,  0.000,   0.500});
            points.push_back({ -0.500,  0.000,   0.000});
            points.push_back({  0.500,  0.000,   0.500});
            points.push_back({ -0.500,  0.000,  -0.500});
            points.push_back({  0.000,  0.500,   0.000});
            points.push_back({1.-zeta, 1.-eta, 1.-zeta});
            points.push_back({   zeta,    eta,    zeta});
            points.push_back({  -zeta, 1.-eta,   -zeta});
            points.push_back({  1.-mu, 1.-eta,  1.-dlt});
            points.push_back({    -mu, 1.-eta,    -dlt});
            points.push_back({ 1.-phi,  0.500,     phi});
            points.push_back({ phi-1.,  0.500,  1.-phi});
            points.push_back({   -dlt, 1.-eta,   1.-mu});
            points.push_back({    dlt,    eta,      mu});
            points.push_back({ psi-1.,  0.000,  1.-psi});
            points.push_back({ 1.-psi,  0.000,     psi});
            points.push_back({   -psi,  0.000,   psi-1});
          }
        }
        if (spgr.get_choice().find("c") != std::string::npos){
          if (2*rlat.get_gamma() < PI){
            if (c_gamma*a/b + s_gamma*s_gamma*a*a/c/c > 1){
              zeta = 0.25/s_gamma/s_gamma + 0.25*a*a/c/c - 0.25*a/b/s_gamma/t_gamma;
              eta = 0.5 + 2.0*b/a*zeta*c_gamma;
              mu = 0.25 + 0.25*a*a/c/c - 0.25/t_gamma/t_gamma + 0.25*b*c_gamma/a*(1.0/s_gamma/s_gamma-a*a/c/c);
              phi = 0.75 - 0.25*a*a/c/c - 0.25*a/t_gamma*(1.0/b/s_gamma-1.0/a/t_gamma);
              psi = 0.5 + (2.0*phi-1.0)*b/a*c_gamma;
              nu = 0.25 + 0.25*(a*a/c/c+a/b/s_gamma/t_gamma-3/t_gamma/t_gamma) + 0.5*b*c_gamma/a*(1/s_gamma/s_gamma-a*a/c/c);
              omega = b*c_gamma/a*(2.0*nu-0.5) + b*s_gamma*t_gamma/a*(2.0*nu-0.5*a*a/c/c-0.5);
              rho = 0.75 + 0.25*c*c/a/a/s_gamma*(a/b/t_gamma-1/s_gamma);
              dlt = (2.0*mu-nu)*b*c_gamma/a + 0.5*omega - 0.25;
              points.push_back({  0.500,  0.500,  0.500});
              points.push_back({  0.500,  0.500,  0.000});
              points.push_back({  0.000,  0.500,  0.000});
              points.push_back({  0.000,  0.000, -0.500});
              points.push_back({  0.000,  0.500, -0.500});
              points.push_back({  0.500,  0.000,  0.000});
              points.push_back({  omega,     nu,     nu});
              points.push_back({1-omega,  1.-nu,  1.-nu});
              points.push_back({  omega,     nu,  nu-1.});
              points.push_back({    eta,   zeta,   zeta});
              points.push_back({ 1.-eta,1.-zeta,  -zeta});
              points.push_back({ 1.-eta,  -zeta,  -zeta});
              points.push_back({  0.500,    rho, 1.-rho});
              points.push_back({  0.500, 1.-rho, rho-1.});
              points.push_back({    psi,    phi, phi-1.});
              points.push_back({ 1.-psi, 1.-phi, 1.-phi});
              points.push_back({ 1.-psi, 1.-phi,   -phi});
              points.push_back({    dlt,     mu,     mu});
              points.push_back({   -dlt,  1.-mu,    -mu});
              points.push_back({   -dlt,    -mu,    -mu});
              points.push_back({    dlt,     mu,   mu-1});
            } else {
              zeta = 0.25/s_gamma/s_gamma + 0.25*a*a/c/c - 0.25*a/b/s_gamma/t_gamma;
              eta = 0.5 + 2*b/a*zeta*c_gamma;
              mu = 0.25 + 0.25*a*a/c/c;
              dlt = (2*mu-0.5)*b/a*c_gamma;
              phi = 0.75 - 0.25*a*a/c/c - 0.25*a/t_gamma*(1/b/s_gamma-1/a/t_gamma);
              psi = 0.5 + (2*phi-1)*b/a*c_gamma;
              points.push_back({  0.500,  0.500, -0.500});
              points.push_back({  0.500,  0.500,  0.000});
              points.push_back({  0.000,  0.500,  0.000});
              points.push_back({  0.000,  0.000, -0.500});
              points.push_back({  0.000,  0.500, -0.500});
              points.push_back({  0.500,  0.000,  0.000});
              points.push_back({ 1.-psi, 1.-phi, 1.-phi});
              points.push_back({    psi,    phi, phi-1.});
              points.push_back({ 1.-psi, 1.-phi,   -phi});
              points.push_back({    eta,   zeta,   zeta});
              points.push_back({ 1.-eta,1.-zeta,  -zeta});
              points.push_back({ 1.-eta,  -zeta,  -zeta});
              points.push_back({    dlt,     mu,     mu});
              points.push_back({   -dlt,  1.-mu,    -mu});
              points.push_back({   -dlt,    -mu,    -mu});
              points.push_back({    dlt,     mu,  mu-1.});
            }
          } else {
            zeta = 0.5/s_gamma/s_gamma - 0.25*a/b/t_gamma/s_gamma;
            eta = 0.5 + 2*b*c_gamma*zeta/a;
            psi = 0.75 - 0.25*c*c/a/a/s_gamma/s_gamma;
            phi = psi + 0.25*c*c/a/b/t_gamma/s_gamma;
            mu = psi + 0.25*(2-a/b/c_gamma)/t_gamma/t_gamma;
            dlt = 1.0 - 0.5*a/t_gamma*(1/b/s_gamma-2/a/t_gamma) - mu;
            points.push_back({  0.000,  0.000,  0.000});
            points.push_back({  0.500,  0.500,  0.500});
            points.push_back({  0.500,  0.500,  0.000});
            points.push_back({  0.000,  0.500,  0.000});
            points.push_back({  0.000,  0.000, -0.500});
            points.push_back({  0.000,  0.500,  0.500});
            points.push_back({  0.000, -0.500, -0.500});
            points.push_back({  0.500,  0.000,  0.000});
            points.push_back({ 1.-eta,1.-zeta,1.-zeta});
            points.push_back({    eta,   zeta,   zeta});
            points.push_back({ 1.-eta,  -zeta,  -zeta});
            points.push_back({ 1.-eta, 1.-dlt,  1.-mu});
            points.push_back({ 1.-eta,   -dlt,    -mu});
            points.push_back({  0.500,    phi, 1.-phi});
            points.push_back({  0.500, 1.-phi, phi-1.});
            points.push_back({ 1.-eta,  1.-mu,   -dlt});
            points.push_back({    eta,     mu,    dlt});
            points.push_back({  0.000, 1.-psi, psi-1.});
            points.push_back({  0.000,    psi, 1.-psi});
            points.push_back({  0.000,  psi-1.,  -psi});
          }
        }
      }
      if (brav == Bravais::A || brav == Bravais::_)
          throw std::logic_error("Unknown monoclinic centring");
    case Holohedry::triclinic:
      // half of the first Brillouin zone is necessary for Triclinic P -1
      // and the whole first Brillouin zone for Triclinic P 1.
    break;
    default:
      throw std::logic_error("Undefined pointgroup holohedry.");
  }
  bool ir_is_ok = this->check_ir_polyhedron();
  // We need to have at least four points to create a polyhedron:
  if (!ir_is_ok && points.size()>3){
    debug_update("Creating an irBZ polyhedron from points.");
    LQVec<double> verts(this->outerlattice, points);
    ir_is_ok = this->set_ir_vertices(verts);
  }
  if (!ir_is_ok){
    LQVec<double> zhat(this->outerlattice, 1u);
    for (size_t i=0; i<3u; ++i) zhat.insert( i==2 ? 1 : 0, 0, i);
    this->set_ir_wedge_normals(zhat);
  }
  return ir_is_ok || this->check_ir_polyhedron();
}
