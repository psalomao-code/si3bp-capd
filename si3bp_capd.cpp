#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "capd/capdlib.h"

using namespace capd;
using capd::autodiff::Node;
using std::cerr;
using std::cout;
using std::endl;

/*
  Rigorous CAPD certificate for the spatial isosceles three-body problem

  Parameters:
      beta = 3/100,
      e    = 987/1000,
      H    = -1.

  Hence
      alpha = 12/97,
      varpi^2 = 25831/1800.

  State ordering: x = (r,z,p_r,p_z).

  The program attempts to certify:

  C1. A symmetric period-four point q of the periapsis map P and
      hyperbolicity of F=P^4.

  C3. A local hyperbolic block and an invariant unstable cone.  The graph
      transform then gives the TRUE local unstable manifold as a graph
      v=h(u), |h'|<=a, without assuming a numerical polynomial ansatz.

  C4. A symmetric homoclinic topological crossing.  Rigorous Poincare
      propagation of the two endpoint tubes proves a sign change of z at the
      next apoapsis.  Subdivided C^1 propagation of the full tangent cone
      proves strict monotonicity of p_z along the unstable arc.  Reversibility
      and the Burns--Weiss criterion then imply positive topological entropy.
      An optional stronger estimate on dz/du certifies transversality (C4+).

  IMPORTANT:
  The energy lift first proves existence AND uniqueness uniformly by endpoint
  signs and H_r<0; interval Newton is used only to narrow the known root.

  This source is written against the public CAPD::DynSys API documented for
  IPoincareMap / C1Rect2Set / C1HORect2Set.  The numerical constants defining
  the small local coordinate change are only *centering data*; every claim is
  accepted only if the interval inequalities printed below pass.
*/

// -----------------------------------------------------------------------------
// Exact problem parameters, represented by outward-rounded CAPD intervals.
// -----------------------------------------------------------------------------

static const interval BETA  = interval(3.) / interval(100.);
static const interval ECC   = interval(987.) / interval(1000.);
static const interval ALPHA = interval(12.) / interval(97.);
static const interval W2    = interval(25831.) / interval(1800.);
static const interval ONE(1.);
static const interval ZERO(0.);

// State-coordinate indices.
enum { R = 0, Z = 1, PR = 2, PZ = 3 };

// -----------------------------------------------------------------------------
// Vector field.
// -----------------------------------------------------------------------------

void si3bpVectorField(Node /*t*/, Node in[], int /*dimIn*/,
                      Node out[], int /*dimOut*/,
                      Node params[], int /*noParams*/)
{
  Node alpha = params[0];
  Node w2    = params[1];
  Node r     = in[R];
  Node z     = in[Z];
  Node pr    = in[PR];
  Node pz    = in[PZ];

  Node c = 1 + 2*alpha;
  Node D = r*r + c*z*z;
  // -3/2 is exactly representable in binary floating point, as recommended
  // in CAPD's own PCR3BP example.
  Node invD32 = D ^ (-1.5);

  out[R]  = pr;
  out[Z]  = pz;
  out[PR] = w2/(r*r*r) - 1/(r*r) - 4*r*invD32/alpha;
  out[PZ] = -4*c*z*invD32/alpha;
}

// -----------------------------------------------------------------------------
// Interval helpers.
// -----------------------------------------------------------------------------

double lo(const interval& x) { return x.leftBound(); }
double hi(const interval& x) { return x.rightBound(); }

double supAbs(const interval& x)
{
  return std::max(std::fabs(lo(x)), std::fabs(hi(x)));
}

double infAbs(const interval& x)
{
  if(lo(x) <= 0.0 && hi(x) >= 0.0) return 0.0;
  return std::min(std::fabs(lo(x)), std::fabs(hi(x)));
}

bool strictlyPositive(const interval& x) { return lo(x) > 0.0; }
bool strictlyNegative(const interval& x) { return hi(x) < 0.0; }
bool containsZero(const interval& x) { return lo(x) <= 0.0 && hi(x) >= 0.0; }

interval intersectIntervals(const interval& a, const interval& b)
{
  double L = std::max(lo(a),lo(b));
  double U = std::min(hi(a),hi(b));
  if(L > U) throw std::runtime_error("empty interval intersection");
  return interval(L,U);
}

interval absInterval(const interval& x)
{
  if(containsZero(x)) return interval(0.0, supAbs(x));
  return interval(infAbs(x), supAbs(x));
}

interval hullInterval(const interval& a, const interval& b)
{
  return interval(std::min(lo(a),lo(b)), std::max(hi(a),hi(b)));
}

// -----------------------------------------------------------------------------
// Energy surface H=-1.
// -----------------------------------------------------------------------------

interval Dpot(const interval& r, const interval& z)
{
  return r*r + (ONE + interval(2.)*ALPHA)*z*z;
}

interval energyPlusOne(const interval& r, const interval& z,
                       const interval& pr, const interval& pz)
{
  interval D = Dpot(r,z);
  return (pr*pr + pz*pz)/interval(2.)
       + W2/(interval(2.)*r*r)
       - ONE/r
       - interval(4.)/(ALPHA*sqrt(D))
       + ONE;
}

interval dHdr(const interval& r, const interval& z)
{
  interval D = Dpot(r,z);
  return -W2/(r*r*r)
       + ONE/(r*r)
       + interval(4.)*r/(ALPHA*D*sqrt(D));
}

interval dHdz(const interval& r, const interval& z)
{
  interval D = Dpot(r,z);
  return interval(4.)*(ONE+interval(2.)*ALPHA)*z
       /(ALPHA*D*sqrt(D));
}

/* Exact explicit lift on z=0, p_r=0, choosing the periapsis branch.

   H=-1 gives
     (w2/2) u^2 - beta^{-1} u + (p_z^2/2+1)=0,  u=1/r.
   The '+' square-root branch gives the smaller r (periapsis).
*/
interval liftRz0(const interval& pz)
{
  interval invBeta = ONE/BETA;
  interval disc = invBeta*invBeta - W2*(pz*pz + interval(2.));
  if(lo(disc) <= 0.0)
    throw std::runtime_error("liftRz0: discriminant is not positive");
  interval u = (invBeta + sqrt(disc))/W2;
  return ONE/u;
}

/* Rigorous parametric lift r=r(p_z,z) on the periapsis energy branch.

   IMPORTANT: interval Newton alone does not prove existence for every
   parameter pair (p_z,z).  We first prove, uniformly on the parameter box,

      Phi(0.36;p_z,z) > 0,   Phi(0.39;p_z,z) < 0,   H_r < 0

   on r in [0.36,0.39].  IVT + strict monotonicity then give exactly one root
   for every parameter pair.  The subsequent interval Newton steps only
   narrow an already known family of roots, hence cannot lose any of them.
*/
interval liftR(const interval& pz, const interval& z)
{
  const interval X0(0.36,0.39);

  interval fL = energyPlusOne(interval(0.36),z,ZERO,pz);
  interval fR = energyPlusOne(interval(0.39),z,ZERO,pz);
  interval d0 = dHdr(X0,z);

  if(!strictlyPositive(fL))
    throw std::runtime_error("liftR: could not certify Phi(0.36)>0 uniformly");
  if(!strictlyNegative(fR))
    throw std::runtime_error("liftR: could not certify Phi(0.39)<0 uniformly");
  if(!strictlyNegative(d0))
    throw std::runtime_error("liftR: could not certify H_r<0 on [0.36,0.39]");

  interval X=X0;
  for(int k=0;k<12;++k){
    double m = 0.5*(lo(X)+hi(X));
    interval M(m);
    interval d = dHdr(X,z);
    if(!strictlyNegative(d))
      throw std::runtime_error("liftR: H_r lost strict negativity during Newton");

    // Parametric interval Newton enclosure.  Since existence and uniqueness
    // were proved above, intersecting with N preserves every exact root
    // r(p_z,z) represented by the input parameter box.
    interval N = M - energyPlusOne(M,z,ZERO,pz)/d;
    interval newX = intersectIntervals(X,N);
    if(hi(newX)-lo(newX) >= hi(X)-lo(X) - 1e-18){
      X = newX;
      break;
    }
    X = newX;
  }

  if(!strictlyNegative(dHdr(X,z)))
    throw std::runtime_error("liftR: final H_r enclosure is not negative");
  return X;
}

IVector liftState(const interval& pz, const interval& z)
{
  IVector x(4);
  x[R]  = (lo(z)==0.0 && hi(z)==0.0) ? liftRz0(pz) : liftR(pz,z);
  x[Z]  = z;
  x[PR] = ZERO;
  x[PZ] = pz;
  return x;
}

// Derivatives of the energy graph r=r(p_z,z).
void liftDerivatives(const interval& r, const interval& pz, const interval& z,
                     interval& rp, interval& rz)
{
  interval hr = dHdr(r,z);
  if(containsZero(hr))
    throw std::runtime_error("liftDerivatives: H_r contains zero");
  rp = -pz/hr;
  rz = -dHdz(r,z)/hr;
}

// -----------------------------------------------------------------------------
// Small matrix helpers.
// -----------------------------------------------------------------------------

IMatrix mul22(const IMatrix& A, const IMatrix& B)
{
  IMatrix C(2,2);
  for(int i=0;i<2;++i)
    for(int j=0;j<2;++j){
      C[i][j]=ZERO;
      for(int k=0;k<2;++k) C[i][j] += A[i][k]*B[k][j];
    }
  return C;
}

IMatrix mul44(const IMatrix& A, const IMatrix& B)
{
  IMatrix C(4,4);
  for(int i=0;i<4;++i)
    for(int j=0;j<4;++j){
      C[i][j]=ZERO;
      for(int k=0;k<4;++k) C[i][j] += A[i][k]*B[k][j];
    }
  return C;
}

IMatrix inverse2(const IMatrix& A)
{
  interval d = A[0][0]*A[1][1]-A[0][1]*A[1][0];
  if(containsZero(d)) throw std::runtime_error("inverse2: singular interval matrix");
  IMatrix B(2,2);
  B[0][0]= A[1][1]/d;
  B[0][1]=-A[0][1]/d;
  B[1][0]=-A[1][0]/d;
  B[1][1]= A[0][0]/d;
  return B;
}

// -----------------------------------------------------------------------------
// Poincare-map evaluation and restriction to the fixed energy surface.
// -----------------------------------------------------------------------------

struct MapEval {
  IVector image;      // full 4D image on p_r=0
  IMatrix D2;         // derivative in section coordinates (p_z,z)
  interval returnTime;

  MapEval() : image(4), D2(2,2), returnTime(0.) {}
};

MapEval evalP(IPoincareMap& pm,
              const interval& pz, const interval& z,
              int iterate)
{
  MapEval out;
  IVector x = liftState(pz,z);

  C1HORect2Set set(x);
  IMatrix mon(4,4);
  interval T;
  IVector y = pm(set,mon,T,iterate);
  IMatrix DPfull = pm.computeDP(y,mon,T);

  interval rp,rz;
  liftDerivatives(x[R],pz,z,rp,rz);

  // B = derivative of lift (p_z,z) -> (r,z,p_r,p_z).
  interval B[4][2];
  for(int i=0;i<4;++i)
    for(int j=0;j<2;++j) B[i][j]=ZERO;
  B[R][0]=rp;   B[R][1]=rz;
  B[Z][0]=ZERO; B[Z][1]=ONE;
  B[PR][0]=ZERO;B[PR][1]=ZERO;
  B[PZ][0]=ONE; B[PZ][1]=ZERO;

  // Output section coordinates are (p_z,z): rows 3 and 1.
  const int row[2] = {PZ,Z};
  for(int i=0;i<2;++i)
    for(int j=0;j<2;++j){
      out.D2[i][j]=ZERO;
      for(int k=0;k<4;++k)
        out.D2[i][j] += DPfull[row[i]][k]*B[k][j];
    }

  out.image = y;
  out.returnTime = T;
  return out;
}


// Hull two MapEval objects.  Used only after a finite subdivision of the
// input rectangle, so it preserves rigor over the union of the cells.
MapEval hullMapEval(const MapEval& A, const MapEval& B)
{
  MapEval H;
  for(int i=0;i<4;++i)
    H.image[i]=hullInterval(A.image[i],B.image[i]);
  for(int i=0;i<2;++i)
    for(int j=0;j<2;++j)
      H.D2[i][j]=hullInterval(A.D2[i][j],B.D2[i][j]);
  H.returnTime=hullInterval(A.returnTime,B.returnTime);
  return H;
}

// One rigorous Poincare return, with adaptive subdivision of the input
// section rectangle if CAPD loses transversality because of wrapping.
// The exact orbit point is contained in one of the cells; hulling the
// validated cells gives a rigorous enclosure over their union.
MapEval evalP1Adaptive(IPoincareMap& pm,
                       const interval& pz, const interval& z,
                       int depth=0)
{
  try{
    return evalP(pm,pz,z,1);
  }catch(const std::exception& e){
    const double pL=lo(pz), pU=hi(pz);
    const double zL=lo(z),  zU=hi(z);
    const double wp=pU-pL, wz=zU-zL;
    if(depth>=24 || (wp<=0.0 && wz<=0.0))
      throw;

    if(depth==0)
      cout << "  one-return enclosure requested subdivision after: "
           << e.what() << "\n";

    if(wp>=wz && wp>0.0){
      const double m=0.5*(pL+pU);
      if(!(pL<m && m<pU)) throw;
      MapEval A=evalP1Adaptive(pm,interval(pL,m),z,depth+1);
      MapEval B=evalP1Adaptive(pm,interval(m,pU),z,depth+1);
      return hullMapEval(A,B);
    }else{
      const double m=0.5*(zL+zU);
      if(!(zL<m && m<zU)) throw;
      MapEval A=evalP1Adaptive(pm,pz,interval(zL,m),depth+1);
      MapEval B=evalP1Adaptive(pm,pz,interval(m,zU),depth+1);
      return hullMapEval(A,B);
    }
  }
}

// Compose several section-to-section returns, but restart from the rigorous
// section-coordinate enclosure after every return.  This removes the large
// artificial time-of-flight wrapping that occurs when asking PoincareMap for
// a high iterate in one shot.
MapEval evalPChainRestart(IPoincareMap& pm,
                          const interval& pz0, const interval& z0,
                          int n)
{
  if(n<1) throw std::runtime_error("evalPChainRestart: n must be positive");

  interval pz=pz0, z=z0;
  IMatrix Dtot(2,2);
  Dtot[0][0]=ONE; Dtot[0][1]=ZERO;
  Dtot[1][0]=ZERO; Dtot[1][1]=ONE;
  interval Ttot(0.);
  MapEval last;

  for(int k=0;k<n;++k){
    MapEval E=evalP1Adaptive(pm,pz,z);
    Dtot=mul22(E.D2,Dtot);
    Ttot += E.returnTime;
    pz=E.image[PZ];
    z =E.image[Z];
    last=E;
  }

  last.D2=Dtot;
  last.returnTime=Ttot;
  return last;
}

/* P^24 followed by the next PlusMinus crossing (apoapsis).
   This is A o F^6 with F=P^4.
*/
MapEval evalToApo(IPoincareMap& pmPeri, IPoincareMap& pmApo,
                  const interval& pz, const interval& z)
{
  // Twenty-four periapsis returns, synchronized one return at a time.
  MapEval P24=evalPChainRestart(pmPeri,pz,z,24);

  // From the final periapsis go to the next apoapsis.  Restarting from the
  // rigorous section-coordinate enclosure remains rigorous because the exact
  // image lies in that rectangle on the fixed energy surface.
  MapEval A=evalP1Adaptive(pmApo,P24.image[PZ],P24.image[Z]);

  MapEval out=A;
  out.D2=mul22(A.D2,P24.D2);
  out.returnTime=P24.returnTime+A.returnTime;
  return out;
}

interval prDot(const IVector& x)
{
  interval D = Dpot(x[R],x[Z]);
  return W2/(x[R]*x[R]*x[R])
       - ONE/(x[R]*x[R])
       - interval(4.)*x[R]/(ALPHA*D*sqrt(D));
}

// -----------------------------------------------------------------------------
// C1: symmetric period-four orbit and hyperbolicity.
// -----------------------------------------------------------------------------

// Value-only evaluation of P.  This deliberately avoids computeDP: for the
// sign bracket and bisection we need only continuity of g(p)=pi_p P(p,0).
// Avoiding computeDP here is much more robust because the derivative formula
// divides by the transverse section speed at the return point.
IVector evalPValue(IPoincareMap& pm,
                   const interval& pz, const interval& z,
                   int iterate)
{
  IVector x = liftState(pz,z);
  C1HORect2Set set(x);
  IMatrix mon(4,4);
  interval T;
  return pm(set,mon,T,iterate);
}

interval gAtPoint(IPoincareMap& pmPeri, double p)
{
  IVector y = evalPValue(pmPeri,interval(p),ZERO,1);
  return y[PZ];
}


// Validate that the first periapsis return is defined and transversal on a
// finite cover of the whole sign bracket.  Successful CAPD Poincare
// evaluations on the cells give the continuity domain needed for the IVT.
bool validateReturnCover(IPoincareMap& pmPeri,
                         const interval& I,
                         int depth=0)
{
  try{
    (void)evalPValue(pmPeri,I,ZERO,1);
    return true;
  }catch(const std::exception& e){
    const double L=lo(I), U=hi(I);
    if(depth>=24 || !(L<U)){
      if(depth==0)
        cout << "  return-cover validation failed: " << e.what() << "\\n";
      return false;
    }
    const double M=0.5*(L+U);
    if(!(L<M && M<U)) return false;
    if(depth==0)
      cout << "  return-cover requested subdivision after: "
           << e.what() << "\\n";
    return validateReturnCover(pmPeri,interval(L,M),depth+1)
        && validateReturnCover(pmPeri,interval(M,U),depth+1);
  }
}

IMatrix hull22(const IMatrix& A, const IMatrix& B)
{
  IMatrix H(2,2);
  for(int i=0;i<2;++i)
    for(int j=0;j<2;++j)
      H[i][j]=hullInterval(A[i][j],B[i][j]);
  return H;
}

// Hyperbolicity of the symmetric period-four point without a four-return
// C^1 integration.  For q=(p_*,0), g(p_*)=0 implies q1=P(q) lies in Fix(S),
// hence P^2(q)=Cq.  Since C(p,z)=(-p,-z), DC=-I and CP=PC give
// DP(Cx)=DP(x) in section coordinates.  If A=DP(q) and B=DP(q1), then
// D(P^4)(q)=B A B A=(B A)^2.
//
// We only need one-return derivative enclosures at q and q1.  The exact q1
// has p_z=0, so after enclosing its z-coordinate we may rigorously restart
// from the rectangle {p_z=0, z in z1}.
IMatrix evalDF4BySymmetry(IPoincareMap& pmPeri,
                          const interval& pRoot,
                          interval& q1Z)
{
  MapEval Aq=evalP1Adaptive(pmPeri,pRoot,ZERO);
  if(!containsZero(Aq.image[PZ]))
    throw std::runtime_error(
      "evalDF4BySymmetry: first-return image does not contain p_z=0");

  q1Z=Aq.image[Z];
  MapEval Bq1=evalP1Adaptive(pmPeri,ZERO,q1Z);

  // Optional consistency checks for P(q1)=Cq at the exact symmetric point.
  interval minusP(-hi(pRoot),-lo(pRoot));
  bool pConsistent = !(hi(Bq1.image[PZ]) < lo(minusP)
                    || lo(Bq1.image[PZ]) > hi(minusP));
  bool zConsistent = containsZero(Bq1.image[Z]);
  if(!(pConsistent && zConsistent))
    throw std::runtime_error(
      "evalDF4BySymmetry: symmetry consistency enclosure failed");

  IMatrix BA=mul22(Bq1.D2,Aq.D2);
  return mul22(BA,BA);
}

struct C1Certificate {
  bool pass=false;
  interval pRoot;
  interval gLeft,gRight;
  IMatrix DF;
  interval lambdaU,lambdaS;

  C1Certificate() : pRoot(0.),gLeft(0.),gRight(0.),DF(2,2),lambdaU(0.),lambdaS(0.) {}
};

C1Certificate certifyC1(IPoincareMap& pmPeri)
{
  C1Certificate C;
  double L=8.5768, U=8.5769;

  C.gLeft  = gAtPoint(pmPeri,L);
  C.gRight = gAtPoint(pmPeri,U);

  cout << "C1 initial bracket:\n";
  cout << "  g(" << std::setprecision(17) << L << ") = " << C.gLeft << "\n";
  cout << "  g(" << U << ") = " << C.gRight << "\n";

  // No monotonicity estimate is needed.  The opposite endpoint signs and
  // continuity already imply at least one p_* in (L,U) with g(p_*)=0.
  if(!(strictlyNegative(C.gLeft) && strictlyPositive(C.gRight))){
    cout << "C1 FAIL: endpoint sign bracket did not validate.\n";
    return C;
  }

  if(!validateReturnCover(pmPeri,interval(L,U))){
    cout << "C1 FAIL: could not validate a continuous transversal first-return "
            "map on the full sign bracket.\n";
    return C;
  }
  cout << "  first return validated on a finite cover of the full sign bracket\n";

  // Rigorous sign bisection.  At every stage the endpoints retain opposite
  // strict signs, so the intermediate value theorem preserves existence of
  // at least one zero inside the current bracket.
  for(int k=0;k<42;++k){
    double m=0.5*(L+U);
    interval gm=gAtPoint(pmPeri,m);
    if(strictlyNegative(gm)) L=m;
    else if(strictlyPositive(gm)) U=m;
    else break; // the enclosure itself already contains the zero
  }
  C.pRoot=interval(L,U);

  // Hyperbolicity of F=P^4 is checked only after the sign bisection has
  // reduced the parameter interval.  This is the first place in C1 where
  // computeDP is genuinely needed.
  cout << "  sign-bisection root bracket width = "
       << (hi(C.pRoot)-lo(C.pRoot)) << "\n";
  try{
    interval q1Z;
    C.DF = evalDF4BySymmetry(pmPeri,C.pRoot,q1Z);
    cout << "  first symmetric image has z in " << q1Z << "\n";
  }catch(const std::exception& e){
    cout << "C1 FAIL during one-return symmetry derivative evaluation: " << e.what() << "\n";
    cout << "  (The existence of the period-four point is already certified; "
            "only the derivative enclosure failed.)\n";
    return C;
  }

  interval tr  = C.DF[0][0]+C.DF[1][1];

  // On the fixed-energy section p_r=0 the restricted symplectic form is
  // dz ^ dp_z (equivalently -dp_z ^ dz in our coordinate order), so every
  // Poincare return has determinant exactly one.  Use this exact identity
  // instead of an interval determinant assembled from correlated entries.
  interval detBox = C.DF[0][0]*C.DF[1][1]-C.DF[0][1]*C.DF[1][0];
  interval disc=tr*tr-interval(4.);
  if(lo(disc)<=0.0){
    cout << "C1 FAIL: trace does not certify hyperbolicity; tr=" << tr
         << " disc=" << disc << "\n";
    return C;
  }
  C.lambdaU=(tr+sqrt(disc))/interval(2.);
  C.lambdaS=(tr-sqrt(disc))/interval(2.);

  cout << "  root p_z in " << C.pRoot << "\n";
  cout << "  DF(q) =\n" << C.DF << "\n";
  cout << "  trace DF = " << tr << "\n";
  cout << "  entrywise determinant enclosure (diagnostic) = " << detBox << "\n";
  cout << "  exact symplectic determinant = 1\n";
  cout << "  lambda_u = " << C.lambdaU << "\n";
  cout << "  lambda_s = " << C.lambdaS << "\n";

  bool eigOK = lo(C.lambdaU)>5.5 && hi(C.lambdaU)<5.8
            && lo(C.lambdaS)>0.17 && hi(C.lambdaS)<0.19;
  C.pass=eigOK;
  cout << (C.pass ? "C1 PASS\n" : "C1 FAIL: eigenvalue target intervals not certified.\n");
  return C;
}

// -----------------------------------------------------------------------------
// Local eigen coordinates and cone/block certificate.
// -----------------------------------------------------------------------------

struct LocalCoords {
  // Columns are approximate unstable/stable eigenvectors.  These decimal
  // constants merely define coordinates; interval tests below absorb error.
  IMatrix B;
  IMatrix Binv;
  LocalCoords() : B(2,2), Binv(2,2)
  {
    B[0][0]=interval( 0.68026252175126);
    B[1][0]=interval(-0.73296855423723);
    B[0][1]=interval( 0.68026252100350);
    B[1][1]=interval( 0.73296855493122);
    Binv=inverse2(B);
  }
};

void localBoxToPhysical(const interval& qP,
                        const LocalCoords& C,
                        const interval& u, const interval& v,
                        interval& p, interval& z)
{
  p = qP + C.B[0][0]*u + C.B[0][1]*v;
  z =       C.B[1][0]*u + C.B[1][1]*v;
}

void physicalToLocal(const interval& qP,
                     const LocalCoords& C,
                     const interval& p, const interval& z,
                     interval& u, interval& v)
{
  interval dp=p-qP;
  u=C.Binv[0][0]*dp+C.Binv[0][1]*z;
  v=C.Binv[1][0]*dp+C.Binv[1][1]*z;
}

IMatrix localJacobian(const IMatrix& D2, const LocalCoords& C)
{
  return mul22(mul22(C.Binv,D2),C.B);
}

struct ConeCertificate {
  bool pass=false;
  double U=0.0;
  double a=0.0;
  IMatrix J;
  interval wholeVImage,rightUImage,leftUImage;
  ConeCertificate() : J(2,2),wholeVImage(0.),rightUImage(0.),leftUImage(0.) {}
};

ConeCertificate tryCone(IPoincareMap& pmPeri,
                        const interval& qP,
                        double U, double a,
                        const LocalCoords& C)
{
  ConeCertificate cert;
  cert.U=U; cert.a=a;
  double V=a*U;

  interval p,z;
  localBoxToPhysical(qP,C,interval(-U,U),interval(-V,V),p,z);
  MapEval E=evalPChainRestart(pmPeri,p,z,4);
  cert.J=localJacobian(E.D2,C);

  // Forward unstable-cone invariance:
  // |J21|+a|J22| <= a( inf|J11| - a sup|J12| ).
  double denom=infAbs(cert.J[0][0]) - a*supAbs(cert.J[0][1]);
  double numer=supAbs(cert.J[1][0]) + a*supAbs(cert.J[1][1]);
  bool cone = denom>1.0 && numer < a*denom;

  // Stable-coordinate trapping for the whole block.
  interval uOut,vOut;
  physicalToLocal(qP,C,E.image[PZ],E.image[Z],uOut,vOut);
  cert.wholeVImage=vOut;
  bool stableTrap=(lo(vOut)>-V && hi(vOut)<V);

  // Exit faces in unstable direction.
  interval pR,zR,pL,zL;
  localBoxToPhysical(qP,C,interval(U),interval(-V,V),pR,zR);
  localBoxToPhysical(qP,C,interval(-U),interval(-V,V),pL,zL);
  MapEval ER=evalPChainRestart(pmPeri,pR,zR,4);
  MapEval EL=evalPChainRestart(pmPeri,pL,zL,4);
  interval ur,vr,ul,vl;
  physicalToLocal(qP,C,ER.image[PZ],ER.image[Z],ur,vr);
  physicalToLocal(qP,C,EL.image[PZ],EL.image[Z],ul,vl);
  cert.rightUImage=ur;
  cert.leftUImage=ul;
  bool exits=(lo(ur)>U && hi(ul)<-U);

  cert.pass=cone && stableTrap && exits;
  return cert;
}

ConeCertificate findCone(IPoincareMap& pmPeri,
                         const interval& qP,
                         const LocalCoords& C)
{
  const std::vector<double> Us={4.0e-5,4.5e-5,5.0e-5,6.0e-5,8.0e-5};
  const std::vector<double> As={1e-4,2e-4,5e-4,1e-3,2e-3,5e-3,1e-2,2e-2,5e-2,1e-1};

  ConeCertificate last;
  for(double U:Us){
    for(double a:As){
      cout << "  trying cone U=" << U << " a=" << a << " ... " << std::flush;
      try{
        ConeCertificate c=tryCone(pmPeri,qP,U,a,C);
        cout << (c.pass ? "PASS" : "no") << "\n";
        if(c.pass){
          cout << "    local J =\n" << c.J << "\n";
          cout << "    v(F(N)) = " << c.wholeVImage << "\n";
          cout << "    u(F(right face)) = " << c.rightUImage << "\n";
          cout << "    u(F(left face))  = " << c.leftUImage << "\n";
          return c;
        }
        last=c;
      }catch(const std::exception& e){
        cout << "exception: " << e.what() << "\n";
      }
    }
  }
  return last;
}

// -----------------------------------------------------------------------------
// C4: sign change and transversality on the true unstable graph.
// -----------------------------------------------------------------------------

struct C4Certificate {
  bool pass=false;             // Burns--Weiss topological crossing certificate
  bool transversePass=false;  // stronger C4+ certificate
  interval zLeft,zRight;
  interval tangentP,tangentZ;
  interval apoPrDot;
  interval apoStateP,apoStateZ,apoStateR;
  C4Certificate() : zLeft(0.),zRight(0.),tangentP(0.),tangentZ(0.),apoPrDot(0.),apoStateP(0.),apoStateZ(0.),apoStateR(0.) {}
};

C4Certificate certifyC4(IPoincareMap& pmPeri, IPoincareMap& pmApo,
                        const interval& qP,
                        const ConeCertificate& cone,
                        const LocalCoords& C)
{
  C4Certificate out;
  if(!cone.pass){
    cout << "C4 SKIPPED: no validated local cone/block.\n";
    return out;
  }

  // Eigen-coordinate values bracketing the exploratory homoclinic crossing.
  const double uL=2.68e-5;
  const double uR=2.70e-5;
  const double a=cone.a;

  auto endpoint=[&](double u)->MapEval{
    interval p,z;
    // Proposition (graph transform) gives |v(u)| <= a|u| on the TRUE W^u.
    localBoxToPhysical(qP,C,interval(u),interval(-a*std::fabs(u),a*std::fabs(u)),p,z);
    return evalToApo(pmPeri,pmApo,p,z);
  };

  MapEval L=endpoint(uL);
  MapEval Rr=endpoint(uR);
  out.zLeft=L.image[Z];
  out.zRight=Rr.image[Z];

  cout << "C4 endpoint images:\n";
  cout << "  z_apo(uL) = " << out.zLeft << "\n";
  cout << "  z_apo(uR) = " << out.zRight << "\n";

  bool signChange = strictlyNegative(out.zLeft) && strictlyPositive(out.zRight);
  if(!signChange){
    cout << "C4 FAIL: endpoint tubes do not rigorously force opposite signs.\n";
    return out;
  }

  /*
     Validate the complete unstable arc by subdivision.  A single interval
     evaluation of A o P^24 on [uL,uR] can suffer severe wrapping.  Every
     cell below still contains the true graph because |v(u)| <= a|u|.

     For each cell, propagate the complete tangent cone (1,m), |m|<=a.
     Burns--Weiss only needs monotonicity of the p_z coordinate, so the
     entropy certificate requires tangentP < -ETA.  tangentZ is accumulated
     only for the optional stronger transverse certificate.
  */
  const int pieces=32;
  const double ETA=50.0;
  bool first=true;
  bool allRegular=true;
  bool allPnegative=true;
  bool allZpositive=true;

  for(int j=0;j<pieces;++j){
    double ua=uL+(uR-uL)*double(j)/double(pieces);
    double ub=uL+(uR-uL)*double(j+1)/double(pieces);
    interval uI(ua,ub);
    // u is positive on this strip, hence |v| <= a*u <= a*ub.
    interval vI(-a*ub,a*ub);

    interval pI,zI;
    localBoxToPhysical(qP,C,uI,vI,pI,zI);
    MapEval E=evalToApo(pmPeri,pmApo,pI,zI);

    interval pdot=prDot(E.image);
    interval m(-a,a);
    interval dpdu=C.B[0][0]+C.B[0][1]*m;
    interval dzdu=C.B[1][0]+C.B[1][1]*m;
    interval tanP=E.D2[0][0]*dpdu+E.D2[0][1]*dzdu;
    interval tanZ=E.D2[1][0]*dpdu+E.D2[1][1]*dzdu;

    if(first){
      out.apoStateP=E.image[PZ];
      out.apoStateZ=E.image[Z];
      out.apoStateR=E.image[R];
      out.apoPrDot=pdot;
      out.tangentP=tanP;
      out.tangentZ=tanZ;
      first=false;
    }else{
      out.apoStateP=hullInterval(out.apoStateP,E.image[PZ]);
      out.apoStateZ=hullInterval(out.apoStateZ,E.image[Z]);
      out.apoStateR=hullInterval(out.apoStateR,E.image[R]);
      out.apoPrDot=hullInterval(out.apoPrDot,pdot);
      out.tangentP=hullInterval(out.tangentP,tanP);
      out.tangentZ=hullInterval(out.tangentZ,tanZ);
    }

    bool reg = hi(pdot)<-0.02;
    bool pneg = hi(tanP)<-ETA;
    bool zpos = lo(tanZ)>1.0e5;
    allRegular = allRegular && reg;
    allPnegative = allPnegative && pneg;
    allZpositive = allZpositive && zpos;

    cout << "  piece " << j+1 << "/" << pieces
         << " u=" << uI
         << "  prdot=" << pdot
         << "  dpz/du=" << tanP
         << "  dz/du=" << tanZ << "\n";
  }

  cout << "  hull apo: r=" << out.apoStateR
       << " z=" << out.apoStateZ << " pz=" << out.apoStateP << "\n";
  cout << "  hull p_r dot = " << out.apoPrDot << "\n";
  cout << "  hull tangent p_z component = " << out.tangentP << "\n";
  cout << "  hull tangent z component   = " << out.tangentZ << "\n";

  out.pass = signChange && allRegular && allPnegative;
  out.transversePass = out.pass && allZpositive;

  if(out.pass){
    cout << "C4 PASS (Burns--Weiss): sign change + regular apoapsis + dp_z/du < -"
         << ETA << ".\n";
    if(out.transversePass)
      cout << "C4+ PASS: dz/du > 1e5 as well; homoclinic is transverse.\n";
    else
      cout << "C4+ not certified (not needed for positive entropy).\n";
  }else{
    cout << "C4 FAIL: topological-crossing certificate not validated.\n";
  }
  return out;
}

// -----------------------------------------------------------------------------
// Main.
// -----------------------------------------------------------------------------

int main()
{
  cout << std::setprecision(17);
  try{
    cout << "SI3BP CAPD rigorous certificate\n";
    cout << "beta=" << BETA << " e=" << ECC
         << " alpha=" << ALPHA << " varpi^2=" << W2 << "\n\n";

    const int dim=4, noParam=2;
    IMap vf(si3bpVectorField,dim,dim,noParam);
    vf.setParameter(0,ALPHA);
    vf.setParameter(1,W2);

    // Separate solvers avoid any accidental step-control interaction between
    // the two requested crossing directions.
    IOdeSolver solverPeri(vf,30);
    IOdeSolver solverApo(vf,30);
    solverPeri.setAbsoluteTolerance(1e-13);
    solverPeri.setRelativeTolerance(1e-13);
    solverApo.setAbsoluteTolerance(1e-13);
    solverApo.setRelativeTolerance(1e-13);

    ICoordinateSection sectionPeri(4,PR);
    ICoordinateSection sectionApo(4,PR);

    // p_r changes minus -> plus at periapsis and plus -> minus at apoapsis.
    IPoincareMap pmPeri(solverPeri,sectionPeri,poincare::MinusPlus);
    IPoincareMap pmApo (solverApo, sectionApo, poincare::PlusMinus);
    pmPeri.setMaxReturnTime(10000.0);
    pmApo.setMaxReturnTime(10000.0);
    pmPeri.setBlowUpMaxNorm(1e7);
    pmApo.setBlowUpMaxNorm(1e7);

    cout << "=== C1 ===\n";
    C1Certificate c1=certifyC1(pmPeri);
    if(!c1.pass){
      cout << "\nOVERALL: FAIL at C1\n";
      return 2;
    }

    cout << "\n=== local unstable block / cone ===\n";
    LocalCoords coords;
    ConeCertificate cone=findCone(pmPeri,c1.pRoot,coords);
    if(!cone.pass){
      cout << "No cone from the built-in search list was validated.\n";
      cout << "This is a numerical-tuning failure, not a mathematical counterexample.\n";
      cout << "Try smaller U, a finer coordinate matrix, or C0/C1 HO set types.\n";
      cout << "\nOVERALL: C1 PASS, C4 NOT CERTIFIED\n";
      return 3;
    }

    cout << "\n=== C4 ===\n";
    C4Certificate c4=certifyC4(pmPeri,pmApo,c1.pRoot,cone,coords);

    cout << "\n=== RESULT ===\n";
    cout << "C1: " << (c1.pass ? "PASS" : "FAIL") << "\n";
    cout << "C4 (Burns--Weiss crossing): " << (c4.pass ? "PASS" : "FAIL") << "\n";
    cout << "C4+ (transverse): " << (c4.transversePass ? "PASS" : "not certified") << "\n";
    if(c1.pass && c4.pass){
      cout << "ALL ENTROPY CERTIFICATES PASSED.\n";
      cout << "The validated graph-transform, reversibility and Burns--Weiss criteria\n"
              "prove a homoclinic topological crossing and positive topological entropy.\n";
      if(c4.transversePass)
        cout << "The stronger transverse-homoclinic certificate also passed.\n";
      return 0;
    }
    return 4;
  }
  catch(const std::exception& e){
    cerr << "CAPD certificate aborted with exception: " << e.what() << "\n";
    return 10;
  }
}
