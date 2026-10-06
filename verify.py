#!/usr/bin/env python3
"""Exact computational checks for the rank-twelve elliptic Calabi--Yau example."""
import argparse
import ast
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import platform
import re
import sys
import tempfile

import sympy as S

VERSION = "1.0.1"
HERE = Path(__file__).resolve().parent
DATA = HERE / "data" / "node_table.json"


def require(condition, message):
    """Checks remain active even when Python is invoked with -O."""
    if not condition:
        raise ValueError(message)


def rational(text):
    require(isinstance(text, str) and re.fullmatch(r"-?\d+(?:/[1-9]\d*)?", text),
            "Expected an exact rational coefficient")
    return S.Rational(text)


def load_rows():
    data = json.loads(DATA.read_text(encoding="utf-8"))
    require(data.get("schema_version") == 1, "Unsupported node-table schema")
    rows = []
    for row in data["rows"]:
        scalars = {}
        for key in ("v", "tau"):
            scalars[key] = rational(row[key]["real"]) + S.I*rational(row[key]["imag"])
        labels = row["sections"]
        require(len(labels) == 4 and all(re.fullmatch(r"P'?([126])", x) for x in labels),
                "Expected four section labels in each row")
        rows.append({**scalars, "sections": labels})
    require(len(rows) == 24, "Expected 24 node-table rows")
    return rows


def latex_scalar(text):
    """Read rational Gaussian numbers; never evaluate arbitrary source text."""
    text = text.strip().strip("$")
    text = re.sub(r"\\frac\{([^{}]*)\}\{([^{}]*)\}", r"((\1)/(\2))", text)
    text = re.sub(r"(\d)\s*i", r"\1*I", text).replace("i", "I")

    def visit(node):
        if isinstance(node, ast.Constant) and type(node.value) is int:
            return S.Integer(node.value)
        if isinstance(node, ast.Name) and node.id == "I":
            return S.I
        if isinstance(node, ast.UnaryOp) and isinstance(node.op, (ast.UAdd, ast.USub)):
            value = visit(node.operand)
            return value if isinstance(node.op, ast.UAdd) else -value
        if isinstance(node, ast.BinOp):
            left, right = visit(node.left), visit(node.right)
            if isinstance(node.op, ast.Add): return left + right
            if isinstance(node.op, ast.Sub): return left - right
            if isinstance(node.op, ast.Mult): return left * right
            if isinstance(node.op, ast.Div) and right != 0: return left / right
        raise ValueError("Unsupported entry in the manuscript node table")

    return S.expand(visit(ast.parse(text, mode="eval").body))


def compare_paper(path, rows):
    source = path.read_text(encoding="utf-8")
    begin, end = "% BEGIN EMBEDDED NODE TABLE", "% END EMBEDDED NODE TABLE"
    require(source.count(begin) == source.count(end) == 1, "Expected one marked node table")
    table = source.split(begin)[1].split(end)[0]
    paper_rows = []
    for line in table.splitlines():
        if not (line.startswith("$") and "&" in line and not line.startswith("$v_0$")):
            continue
        cells = line.rstrip(" \\").split("&")
        require(len(cells) == 6, "Expected six cells in the manuscript table")
        labels = []
        for cell in cells[2:]:
            number = re.search(r"\{([126])\}", cell)
            require(number is not None, "Unknown section in the manuscript table")
            labels.append(("P'" if "\\Pp" in cell else "P") + number.group(1))
        paper_rows.append({"v": latex_scalar(cells[0]), "tau": latex_scalar(cells[1]),
                           "sections": labels})
    require(paper_rows == rows, "The manuscript node table differs from data/node_table.json")
    return hashlib.sha256(path.read_bytes()).hexdigest()


def compute(rows):
    v, t, u, z = S.symbols('v t u z')
    i = S.I
    A = lambda x: x**8 + 14*x**4 + 1
    B = lambda x: x*(x**4-1)
    C = lambda x: x**12-33*x**8-33*x**4+1
    a, b, c = A(2), B(2), C(2)
    checks = []

    def zero(expr, name):
        require(S.expand(expr) == 0, name)
        checks.append(name)

    def coeffs(u):
        aa = [-18*(1+i)*(u*v+i),18*(1+i)*(u*v+i),
              18*(1-i)*(u*v-i),18*(1-i)*(u+i*v),-18*(1+i)*(u-i*v),S.Integer(0)]
        bb = [-3*(5*u**4*v**4-u**4-v**4+24*i*u**3*v**3-24*u**2*v**2-24*i*u*v+5)]*2
        bb += [-3*(5*u**4*v**4-u**4-v**4-24*i*u**3*v**3-24*u**2*v**2+24*i*u*v+5),
               3*(u**4*v**4-5*u**4-5*v**4-24*i*u**3*v+24*u**2*v**2+24*i*u*v**3+1),
               3*(u**4*v**4-5*u**4-5*v**4+24*i*u**3*v+24*u**2*v**2-24*i*u*v**3+1),
               3*(u**4*v**4+u**4+v**4+6*u**4*v**2+6*u**2*v**4-12*u**2*v**2+6*u**2+6*v**2+1)]
        return [(aj,bj,S.expand(aj**3/108),S.expand((3*aj**2*bj-11664*B(u)*B(v)-(aj**3/108)**2)/108)) for aj,bj in zip(aa,bb)]

    zero(A(v)**3-C(v)**2-108*B(v)**4, 'ABC identity')
    for j,(aj,bj,cj,dj) in enumerate(coeffs(u),1):
        m=B(u)*B(v)
        zero((54*z**2+cj*z+dj)**2-((aj*z+bj)**3-27*A(u)*A(v)*(aj*z+bj)+2916*(z**4-4*m*z**2+2*m**2)+54*C(u)*C(v)),f'generic section identity {j}')
    points=[]
    for uv in [S.Integer(2),2*i]:
        for aj,bj,cj,dj in coeffs(uv):
            zz=B(uv)*t+B(v)/t
            points.append((S.expand(t**2*(aj*zz+bj)),S.expand(t**3*(54*zz**2+cj*zz+dj))))
    f=-27*a*A(v)*t**4
    g=2916*(b**4*t**10+B(v)**4*t**2)+54*c*C(v)*t**6
    for j,(x,y) in enumerate(points):
        zero(y*y-x**3-f*x-g,f'Weierstrass section {j+1}')
        zero(S.expand(y).coeff(t,1)-54*B(v)**2,f'section {j+1} component at zero')
        zero(S.expand(y).coeff(t,5)-54*b**2*(1 if j<6 else -1),f'section {j+1} component at infinity')
    print('Section identities passed',flush=True)

    J=S.expand(a**3*B(v)**4-b**4*A(v)**3)
    zero(J+(v**4-16)*(16*v**4-1)*(v**4-81)*(81*v**4-1)*(25*v**4-14*v**2+25)*(25*v**4+14*v**2+25),'J factorization')
    require(S.degree(S.gcd(J,S.diff(J,v)),v)==0, "J has no repeated roots")
    for pol in [A(v),B(v),C(v)]:require(S.degree(S.gcd(J,pol),v)==0, "J coprime to A, B and C")
    for eps in [1,-1]:
        zero((c*C(v)+108*eps*b*b*B(v)**2)**2-(a*A(v))**3+108*(b*b*C(v)-eps*c*B(v)**2)**2,f'node identity {eps}')
    zero(J+(b*b*C(v)-c*B(v)**2)*(b*b*C(v)+c*B(v)**2),'J=L product')
    zero(S.rem((c*C(v))**2-4*54*b**4*54*B(v)**4+108*a**3*B(v)**4,A(v),v),'cusp discriminant')

    # Projection quartic and its invariants.
    alpha=(v*v+1)**2; beta=(v*v-1)**2
    k=S.symbols('k')
    pq=S.Poly(((alpha+beta)*k*k-34*t**4)**2-4*(k**3-t**4)*(alpha*beta*k-225*t**4),k)
    q4,q3,q2,q1,q0=pq.all_coeffs()
    II=12*q4*q0-3*q3*q1+q2*q2
    JJ=72*q4*q2*q0+9*q3*q2*q1-27*q4*q1*q1-27*q3*q3*q0-2*q2**3
    zero(II-16*a*A(v)*t**8,'quartic invariant I')
    zero(JJ+128*t**12*(54*b**4*t**4+c*C(v)+54*B(v)**4/t**4),'quartic invariant J')
    w=S.symbols('w')
    zero(w**20*f.subs({v:1/w,t:t/w**3},simultaneous=True)-f.subs(v,w),'f transition')
    zero(w**30*g.subs({v:1/w,t:t/w**3},simultaneous=True)-g.subs(v,w),'g transition')

    # Verify the general pointed quartic-to-Weierstrass transformation.
    ss, eta, qq, a4, a3, a2, a1 = S.symbols('ss eta qq a4 a3 a2 a1')
    uu=(2*qq*(eta+qq)+a1*ss)/ss**2
    zz=((uu**2-4*qq**2*a4)*ss-(a1*uu+2*qq**2*a3))/(2*qq)
    rhs=uu**3+a2*uu**2+(a3*a1-4*a4*qq**2)*uu+a4*a1**2+qq**2*a3**2-4*a4*a2*qq**2
    numerator=S.together(zz**2-rhs).as_numer_denom()[0]
    relation=eta**2-(a4*ss**4+a3*ss**3+a2*ss**2+a1*ss+qq**2)
    zero(S.rem(numerator,relation,eta),'pointed quartic birational transformation')
    invI=12*a4*qq**2-3*a3*a1+a2**2
    invJ=72*a4*a2*qq**2+9*a3*a2*a1-27*a4*a1**2-27*a3**2*qq**2-2*a2**3
    uu0=S.symbols('uu0')
    rhs0=uu0**3+a2*uu0**2+(a3*a1-4*a4*qq**2)*uu0+a4*a1**2+qq**2*a3**2-4*a4*a2*qq**2
    zero(729*rhs0-((9*uu0+3*a2)**3-27*invI*(9*uu0+3*a2)-27*invJ),'pointed quartic invariant scaling')

    # Lowest weighted forms at the six quartic singularities.
    xx, yy, ss, zz = S.symbols('xx yy ss zz')
    binary=xx*yy*(xx+25*yy)*(xx+9*yy)
    def initial(expr):
        p=S.Poly(S.expand(expr),xx,yy,ss,zz)
        weights=[sum(a*b for a,b in zip(mon,[1,1,1,2])) for mon,co in p.terms()]
        low=min(weights)
        require(low==4, "lowest weighted degree is four")
        return sum(co*xx**mon[0]*yy**mon[1]*ss**mon[2]*zz**mon[3] for (mon,co),weight in zip(p.terms(),weights) if weight==low)
    zero(initial((zz-1-ss**4)*(zz**2-4*ss**4)-binary)-(-zz**2+4*ss**4-binary),'weighted form at zero and infinity')
    for vv in [1,-1,i,-i]:
        al=((vv+ss)**2+1)**2;be=((vv+ss)**2-1)**2
        if vv in [i,-i]:al,be=be,al
        zero(initial((zz+al-be/2)*(zz**2-be**2/4)-binary)-(4*zz**2-16*ss**4-binary),f'weighted form at {vv}')

    def order(poly):
        return min(m[0] for m,co in poly.terms() if co)

    M=S.Matrix([[8,-4,2,2,2,2],[-4,8,2,2,2,2],[2,2,8,2,2,2],[2,2,2,8,2,2],[2,2,2,2,8,2],[2,2,2,2,2,8]])
    N=S.Matrix([[3,-3,-3,0,0,0],[-3,3,3,0,0,0],[3,-3,3,0,0,0],[0,0,0,3,3,0],[0,0,0,-3,3,0],[0,0,0,0,0,0]])
    expected=M.row_join(N).col_join(N.T.row_join(M))/3
    height_results={}
    for vv in [5,7]:
        ps=[(S.Poly(x.subs(v,vv),t,extension=i),S.Poly(y.subs(v,vv),t,extension=i)) for x,y in points]
        H=S.eye(12)*S.Rational(8,3)
        intersections=S.zeros(12)
        boundary=[]
        for j in range(12):
            for l in range(j):
                dx,dy=ps[j][0]-ps[l][0],ps[j][1]-ps[l][1]
                gcd=S.gcd(dx,dy)
                interior=gcd.degree()-order(gcd)
                at0=order(dx)-1
                atinf=(3-dx.degree()) if (j<6)==(l<6) else 0
                require(at0>=0 and atinf>=0, "nonnegative boundary intersections")
                inter=interior+at0+atinf
                intersections[j,l]=intersections[l,j]=inter
                correction=S.Rational(2,3)+(S.Rational(2,3) if (j<6)==(l<6) else S.Rational(1,3))
                H[j,l]=H[l,j]=2-inter-correction
                if at0 or atinf:boundary.append([l+1,j+1,at0,atinf])
        require(H==expected, f"height matrix at {vv}: {H}")
        minors=[H[:j,:j].det() for j in range(1,13)]
        require(all(n>0 for n in minors), "positive leading principal minors")
        require(H.det()==S.Rational(4096,9), "height determinant")
        height_results[str(vv)]={'matrix':[[str(x) for x in row] for row in H.tolist()], 'intersection_matrix':intersections.tolist(),'boundary_intersections':boundary,'leading_minors':list(map(str,minors)),'determinant':str(H.det())}
        print(f'Height matrix at v={vv} passed',flush=True)

    P=[(S.Poly(x,v,t,extension=i),S.Poly(y,v,t,extension=i)) for x,y in points]
    fp=S.Poly(f,v,t,extension=i)
    gp=S.Poly(g,v,t,extension=i)
    node_data=[]
    roots=[]
    def ev(poly,vv,tt):return S.expand(poly.eval(v,vv).eval(t,tt))
    for row in rows:
        vv,tau=row["v"],row["tau"];roots.append(vv)
        zero(J.subs(v,vv),f'J root {vv}')
        zero(tau**4-C(vv)/c,f'fourth root {vv}')
        for jj,label in enumerate(row["sections"]):
            number=int(label[-1])
            idx=number-1+(6 if "\'" in label else 0)
            tt=S.expand(i**jj*tau)
            xx=S.cancel(3*a*a*C(vv)/(c*A(vv))*tt**2).expand(complex=True)
            zero(ev(P[idx][0],vv,tt)-xx,f'node x {vv} {jj}')
            zero(ev(P[idx][1],vv,tt),f'node y {vv} {jj}')
            zero(3*xx**2+ev(fp,vv,tt),f'node fiber derivative {vv} {jj}')
            # Independence of dq and d(3p^2+f) proves the ODP condition
            # and smoothness of the section blow-up at each covered node.
            px,py=P[idx]
            h_v=6*xx*ev(px.diff(v),vv,tt)+ev(fp.diff(v),vv,tt)
            h_t=6*xx*ev(px.diff(t),vv,tt)+ev(fp.diff(t),vv,tt)
            det=S.expand(ev(py.diff(v),vv,tt)*h_t-ev(py.diff(t),vv,tt)*h_v)
            require(det!=0, f"nonzero graph Jacobian at {vv}, {tt}, section {idx+1}")
            node_data.append({'v':str(vv),'t':str(tt),'section':("P'" if idx>=6 else 'P')+str(idx%6+1),'graph_jacobian':str(det)})
    require(len(set(roots))==24 and len(node_data)==96, "24 distinct roots and 96 nodes")
    print('All 96 incidence and graph Jacobian checks passed',flush=True)

    # Intersection form in the basis (K_B,C_+,C_-).
    Q=S.Matrix([[2,1,1],[1,-3,0],[1,0,-3]])
    D=S.Matrix([-12,-4,-4]); KK=S.Matrix([1,0,0])
    dsq=(D.T*Q*D)[0]; kd=(KK.T*Q*D)[0]
    pa=1+(dsq+kd)/2
    eD=2-2*pa+96+2*64
    eX=16+(eD-96-64)+192+128
    require((dsq,kd,pa,eD,eX)==(384,-32,177,-128,48), "intersection and Euler arithmetic")
    return {
        "status": "PASS",
        "identity_count": len(checks),
        "identity_checks": checks,
        "heights": height_results,
        "nodes": node_data,
        "invariants": {"D_squared": dsq, "K_D": kd, "arithmetic_genus_D": pa,
                       "Euler_D": eD, "Euler_X": eX},
        "geometric_inputs": {"intersection_matrix": Q.tolist(),
                             "D_coordinates": list(D), "nodes": 96, "cusps": 64,
                             "boundary_euler_contribution": 16},
        "hodge_numbers": {"h11_from_manuscript": 25,
                          "h21_from_h11_and_Euler_X": 25-eX/2},
        "scope": ("Exact algebraic identities, specialized height/intersection calculations, "
                  "and node-incidence and graph-Jacobian checks. Euler arithmetic uses the "
                  "geometric stratification and intersection form proved in the manuscript. "
                  "h11=25 is a manuscript input; h21 is derived from it and Euler_X. "
                  "The written geometric arguments are not formally verified by this program."),
    }


def write_report(path, report):
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8", dir=path.parent,
                                     prefix=".mw-rank12-", delete=False) as handle:
        temporary = Path(handle.name)
        json.dump(report, handle, indent=2, default=str)
        handle.write("\n")
    temporary.replace(path)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path,
                        default=HERE.parent / "MW-rank12-results" / "results.json",
                        help="JSON report path outside the repository; must not overwrite "
                             "--paper (default: %(default)s)")
    parser.add_argument("--paper", type=Path,
                        help="optional manuscript .tex file whose marked node table must match")
    args = parser.parse_args()
    requested_output = Path(os.path.abspath(args.output))
    args.output = requested_output.resolve()
    if any(path == HERE or HERE in path.parents for path in (requested_output, args.output)):
        parser.error("the output must be outside the repository and must not overwrite its files")
    if args.paper and args.output == args.paper.resolve():
        parser.error("the output must not overwrite the manuscript")
    metadata = {
        "software": "MW-rank12", "version": VERSION,
        "python_version": platform.python_version(), "sympy_version": S.__version__,
        "checked_at_utc": datetime.now(timezone.utc).isoformat(),
        "paper_table_compared": args.paper is not None,
    }
    # Invalidate an earlier PASS before work that may fail or be interrupted.
    write_report(args.output, {**metadata, "status": "RUNNING"})
    try:
        metadata["verifier_sha256"] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
        metadata["data_sha256"] = hashlib.sha256(DATA.read_bytes()).hexdigest()
        rows = load_rows()
        if args.paper:
            metadata["source_sha256"] = compare_paper(args.paper, rows)
        report = {**metadata, **compute(rows)}
    except KeyboardInterrupt:
        write_report(args.output, {**metadata, "status": "INTERRUPTED",
                                   "error": "Verification interrupted before completion"})
        print("INTERRUPTED: verification did not complete", file=sys.stderr, flush=True)
        return 130
    except Exception as error:
        write_report(args.output, {**metadata, "status": "FAIL", "error": str(error)})
        print(f"FAIL: {error}", file=sys.stderr, flush=True)
        return 1
    write_report(args.output, report)
    print(f"PASS: {report['identity_count']} identities, two height matrices, 96 node "
          "incidences and nonzero graph Jacobians; e(X)=48.", flush=True)
    print(f"Report: {args.output}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
