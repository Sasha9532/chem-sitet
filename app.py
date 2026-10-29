# -*- coding: utf-8 -*-
# Химический конструктор реакций (10-11 класс). Запуск: python app.py -> http://127.0.0.1:5000
import re
import sys
import subprocess
from math import gcd

try:
    from flask import Flask, request, jsonify, Response
except ImportError:  # автоустановка Flask, если его нет
    subprocess.check_call([sys.executable, "-m", "pip", "install", "flask"])
    from flask import Flask, request, jsonify, Response

app = Flask(__name__)

ROOT = ["", "мет", "эт", "проп", "бут", "пент", "гекс"]
ALK = ["", "метил", "этил", "пропил", "бутил", "пентил", "гексил"]
ACID_T = ["", "муравьиная", "уксусная", "пропионовая", "масляная", "валериановая", "капроновая"]
SALT_T = ["", "формиат", "ацетат", "пропионат", "бутират", "валерат", "капронат"]
ALD_T = ["", "формальдегид", "ацетальдегид", "пропионовый альдегид", "масляный альдегид", "", ""]
HALN = {"Cl": "хлор", "Br": "бром"}
SUBT = str.maketrans("₀₁₂₃₄₅₆₇₈₉−–—", "0123456789---")
LOOK = str.maketrans("СНОКВ", "CHOKB")


def nk(s):
    return re.sub(r"[\s\-()]", "", s)


def qn(s):
    s = re.sub(r"\s+", " ", s.strip().lower().replace("ё", "е"))
    s = re.sub(r"\s*-\s*", "-", s)
    return re.sub(r"(?<=[а-я])\s(?=\d)", "-", s)


def alkyl(n):
    return "CH3" if n == 1 else "C%dH%d" % (n, 2 * n + 1)


# ---------- строение и формулы органических веществ ----------
def spec(c, n, p):
    m = n - 1 if c in ("aldehyde", "acid") else n
    b = [1] * max(m - 1, 0)
    s = [[] for _ in range(m)]
    end = None
    if c == "alkene": b[p - 1] = 2
    elif c == "alkyne": b[p - 1] = 3
    elif c == "alcohol": s[p - 1] = ["OH"]
    elif c == "diol": s[p - 1] = ["OH"]; s[p] = ["OH"]
    elif c == "ketone": s[p - 1] = ["=O"]
    elif c.startswith("dihal"): s[p - 1] = [c[5:]]; s[p] = [c[5:]]
    elif c.startswith("hal"): s[p - 1] = [c[3:]]
    elif c == "aldehyde": end = "CHO"
    elif c == "acid": end = "COOH"
    return m, b, s, end


def render(c, n, p, rev=False, att=True):
    if c == "alkoxide":
        return alkyl(n) + "-ONa"
    if c in ("salt", "ester"):
        a = render("acid", n, 1, rev, att)
        if c == "salt":
            return a.replace("COOH", "COONa").replace("HOOC", "NaOOC")
        return a.replace("COOH", "COO-" + alkyl(p)).replace("HOOC", alkyl(p) + "-OOC")
    m, b, s, end = spec(c, n, p)
    if m == 0:
        return "HCHO" if c == "aldehyde" else "HCOOH"
    A = []
    for i in range(m):
        l = b[i - 1] if i else 0
        r = b[i] if i < m - 1 else (1 if end else 0)
        sb = s[i]
        if "=O" in sb:
            t = "CO"; sb = [x for x in sb if x != "=O"]
        else:
            h = 4 - l - r - len(sb)
            t = "C" + ("H" + (str(h) if h > 1 else "") if h else "")
        A.append([t, sb])
    bonds = list(b)
    if end:
        A.append([end, []]); bonds.append(1)
    if rev:
        A.reverse(); bonds.reverse()
        if end: A[0][0] = {"CHO": "OHC", "COOH": "HOOC"}[end]
    L = len(A); sym = {1: "-", 2: "=", 3: "≡"}; out = ""
    for i, (t, sb) in enumerate(A):
        for x in sb:
            if L == 1: t = t + x if att else t + "-" + x
            elif 0 < i < L - 1: t += "(" + x + ")"
            elif i == 0: t = t + x if att else ("HO" if x == "OH" else x) + "-" + t
            else: t = t + x if att else t + "-" + x
        out += t + (sym[bonds[i]] if i < L - 1 else "")
    return out


def names(c, n, p):
    R, a, t, al = ROOT[n], ALK[n], "", []
    if c == "alkane": i = R + "ан"
    elif c == "alkene":
        i = R + "ен" + ("-%d" % p if n > 3 else ""); t = {2: "этилен", 3: "пропилен"}.get(n, ""); al = ["%s-%d-ен" % (R, p)]
    elif c == "alkyne":
        i = R + "ин" + ("-%d" % p if n > 3 else ""); t = "ацетилен" if n == 2 else ""
    elif c == "alcohol":
        i = R + "анол" + ("-%d" % p if n > 2 else "")
        t = {(1, 1): "метиловый спирт", (2, 1): "этиловый спирт", (3, 1): "пропиловый спирт",
             (3, 2): "изопропиловый спирт", (4, 1): "бутиловый спирт"}.get((n, p), "")
        al = ["%sан-%d-ол" % (R, p), "%d-%sанол" % (p, R)]
    elif c == "diol":
        i = "%sандиол-%d,%d" % (R, p, p + 1); t = {2: "этиленгликоль", 3: "пропиленгликоль"}.get(n, "") if p == 1 else ""
    elif c == "aldehyde":
        i = R + "аналь"; t = ALD_T[n]; al = {1: ["муравьиный альдегид"], 2: ["уксусный альдегид"]}.get(n, [])
    elif c == "ketone":
        i = R + "анон" + ("-%d" % p if n > 4 else ""); t = "ацетон" if n == 3 else ""
    elif c == "acid":
        i = R + "ановая кислота"; t = ACID_T[n] + " кислота"
    elif c.startswith("dihal"):
        i = "%d,%d-ди%s%sан" % (p, p + 1, HALN[c[5:]], R)
    elif c.startswith("hal"):
        i = ("%d-" % p if n > 2 else "") + HALN[c[3:]] + R + "ан"
    elif c == "salt":
        i = R + "аноат натрия"; t = SALT_T[n] + " натрия"
    elif c == "alkoxide":
        i = R + "анолят натрия"; t = a + "ат натрия"
    else:  # ester
        i = ALK[p] + ROOT[n] + "аноат"; t = ALK[p] + SALT_T[n]
    return i, t, al


INF, NAME, FORM = {}, {}, {}
ALL = []


def add(c, n, p=1):
    ALL.append(("o", c, n, p))


for n in range(1, 7):
    for c in ("alkane", "aldehyde", "acid", "salt", "alkoxide"): add(c, n)
    for X in HALN:
        for p in range(1, (n + 1) // 2 + 1): add("hal" + X, n, p)
        if n > 1:
            for p in range(1, n // 2 + 1): add("dihal" + X, n, p)
    if n > 1:
        for p in range(1, n // 2 + 1): add("alkene", n, p); add("alkyne", n, p); add("diol", n, p)
    for p in range(1, (n + 1) // 2 + 1): add("alcohol", n, p)
    if n > 2:
        for p in range(2, (n + 1) // 2 + 1): add("ketone", n, p)
for x in range(1, 5):
    for y in range(1, 5): add("ester", x, y)

# неорганика: ключ -> "название|синонимы"
INORG = {
    "H2": "водород", "O2": "кислород", "Cl2": "хлор", "Br2": "бром|бромная вода", "HCl": "хлороводород (соляная кислота)|соляная кислота",
    "HBr": "бромоводород|бромоводородная кислота", "H2O": "вода", "H2SO4": "серная кислота", "HNO3": "азотная кислота",
    "NaOH": "гидроксид натрия (едкий натр)|едкий натр|щелочь", "KOH": "гидроксид калия (едкое кали)|едкое кали",
    "Na": "натрий", "Cu(OH)2": "гидроксид меди(II)", "CuO": "оксид меди(II)|оксид меди", "Cu2O": "оксид меди(I)", "Cu": "медь",
    "Ag": "серебро", "Ag2O": "оксид серебра (аммиачный раствор)|аммиачный раствор оксида серебра|реактив толленса",
    "KMnO4": "перманганат калия|марганцовка", "MnO2": "оксид марганца(IV)", "Na2CO3": "карбонат натрия|сода",
    "NaHCO3": "гидрокарбонат натрия|питьевая сода", "CO2": "углекислый газ|оксид углерода(IV)", "NH3": "аммиак",
    "NH4Cl": "хлорид аммония|нашатырь", "CaC2": "карбид кальция", "Ca(OH)2": "гидроксид кальция (гашёная известь)|гашеная известь",
    "NaCl": "хлорид натрия (поваренная соль)|поваренная соль", "KCl": "хлорид калия", "NaBr": "бромид натрия", "KBr": "бромид калия",
    "Na2SO4": "сульфат натрия", "NaNO3": "нитрат натрия", "CuCl2": "хлорид меди(II)", "C6H6": "бензол",
    "CH2=CH-Cl": "хлорэтен (винилхлорид)|винилхлорид", "(-CH2-CH2-)n": "полиэтилен", "(-CH2-CH(CH3)-)n": "полипропилен",
}

for s in ALL:
    _, c, n, p = s
    disp = render(c, n, p, c in ("alcohol", "halCl", "halBr") and p == 1, c.startswith(("hal", "dihal")))
    i, t, al = names(c, n, p)
    INF[s] = (disp, i + (" (" + t + ")" if t else ""))
    for rev in (False, True):
        for att in (False, True):
            FORM.setdefault(nk(render(c, n, p, rev, att)), s)
    mo = None
    if p == 1:
        mo = {"alkane": "C%dH%d" % (n, 2 * n + 2), "alkene": "C%dH%d" % (n, 2 * n),
              "alkyne": "C%dH%d" % (n, 2 * n - 2), "alcohol": "C%dH%dOH" % (n, 2 * n + 1)}.get(c)
        if c == "acid" and n > 2: mo = "C%dH%dCOOH" % (n - 1, 2 * n - 1)
    if mo: FORM.setdefault(mo, s)
    for nm in [i, t] + al:
        if nm: NAME.setdefault(qn(nm), s)
    if p == 1:
        NAME.setdefault(qn(re.sub(r"-\d+(,\d+)?$", "", i)), s)
        NAME.setdefault(qn(re.sub(r"^\d+(,\d+)?-", "", i)), s)

for key, val in INORG.items():
    s = ("i", key)
    parts = val.split("|")
    INF[s] = (key, parts[0])
    FORM[nk(key)] = s
    for part in parts:
        NAME.setdefault(qn(part), s)
        NAME.setdefault(qn(re.sub(r"\s+\(.*", "", part)), s)

FORMLOW = {}
for k, v in FORM.items():
    FORMLOW[k.lower()] = None if k.lower() in FORMLOW and FORMLOW[k.lower()] != v else v


def I(k): return ("i", k)


def S(c, n, p=1):
    if c in ("alkene", "alkyne", "diol") or c.startswith("dihal"): p = min(p, n - p)
    elif c in ("alcohol", "ketone") or (c.startswith("hal")): p = min(p, n + 1 - p)
    return ("o", c, n, p)


def mf(a):
    _, c, n, p = a
    return {"alkane": (n, 2 * n + 2, 0), "alkene": (n, 2 * n, 0), "alkyne": (n, 2 * n - 2, 0), "alcohol": (n, 2 * n + 2, 1),
            "diol": (n, 2 * n + 2, 2), "aldehyde": (n, 2 * n, 1), "ketone": (n, 2 * n, 1), "acid": (n, 2 * n, 2),
            "ester": (n + p, 2 * (n + p), 2)}.get(c)


def mk(t, cond, L, R):
    if any(s not in INF for _, s in L + R): return None
    f = lambda X: " + ".join((c if isinstance(c, str) else ("" if c == 1 else str(c))) + INF[s][0] for c, s in X)
    seen, pr = [], []
    for _, s in R:
        if s not in seen:
            seen.append(s); pr.append({"f": INF[s][0], "n": INF[s][1]})
    return {"type": t, "cond": cond, "eq": f(L) + " → " + f(R), "prods": pr}


# ---------- неорганические реакции ----------
RX = """HCl+NaOH>NaCl+H2O|нейтрализация
HCl+KOH>KCl+H2O|нейтрализация
HBr+NaOH>NaBr+H2O|нейтрализация
H2SO4+2NaOH>Na2SO4+2H2O|нейтрализация
HNO3+NaOH>NaNO3+H2O|нейтрализация
Na2CO3+2HCl>2NaCl+CO2+H2O|обмен
NaHCO3+HCl>NaCl+CO2+H2O|обмен
2Na+2H2O>2NaOH+H2|замещение
2Na+2HCl>2NaCl+H2|замещение
H2+Cl2>2HCl|соединение|hν (свет)
H2+Br2>2HBr|соединение|t
CaC2+2H2O>C2H2+Ca(OH)2|гидролиз карбида (получение ацетилена)
NH3+HCl>NH4Cl|соединение
CuO+H2>Cu+H2O|восстановление|t
CuO+2HCl>CuCl2+H2O|обмен
Cu(OH)2+2HCl>CuCl2+2H2O|обмен
CO2+2NaOH>Na2CO3+H2O|обмен""".split("\n")


def lookup(tok):
    m = re.match(r"^(\d*)(.+)$", tok)
    return (int(m.group(1) or 1), FORM[nk(m.group(2))])


RXD = {}
for line in RX:
    eq, typ, *cd = line.split("|")
    l, r = eq.split(">")
    L = [lookup(x) for x in l.split("+")]
    R = [lookup(x) for x in r.split("+")]
    RXD[frozenset(s[1] for _, s in L)] = (typ, cd[0] if cd else "", L, R)


# ---------- органические реакции ----------
def org(a, b):
    _, c, n, p = a
    R = []
    k = b[1] if b and b[0] == "i" else None
    X = k[:-1] if k in ("Cl2", "Br2") else (k[1:] if k in ("HCl", "HBr") else None)

    def ad(*x):
        r = mk(*x)
        if r: R.append(r)

    if k == "O2" and mf(a):
        x, y, z = mf(a)
        o = 4 * x + y - 2 * z
        g = gcd(gcd(4, o), gcd(4 * x, 2 * y))
        ad("горение (полное окисление)", "t", [(4 // g, a), (o // g, b)], [(4 * x // g, I("CO2")), (2 * y // g, I("H2O"))])
    if c == "alkane":
        if k in ("Cl2", "Br2"):
            ad("замещение (радикальное)", "hν (свет)", [(1, a), (1, b)], [(1, S("hal" + X, n, 1 if n < 3 else 2)), (1, I("H" + X))])
        if b is None:
            if n == 1: ad("пиролиз метана", "1500 °C", [(2, a)], [(1, S("alkyne", 2)), (3, I("H2"))])
            else: ad("дегидрирование", "Ni (или Cr2O3), t", [(1, a)], [(1, S("alkene", n)), (1, I("H2"))])
    if c == "alkene":
        q = 2 if (p == 1 and n > 2) else p
        if k == "H2": ad("присоединение (гидрирование)", "Ni, t", [(1, a), (1, b)], [(1, S("alkane", n))])
        if k in ("Cl2", "Br2"): ad("присоединение (галогенирование)", "обесцвечивает бромную воду" if X == "Br" else "CCl4", [(1, a), (1, b)], [(1, S("dihal" + X, n, p))])
        if k in ("HCl", "HBr"): ad("присоединение (гидрогалогенирование, правило Марковникова)", "", [(1, a), (1, b)], [(1, S("hal" + X, n, q))])
        if k == "H2O": ad("присоединение (гидратация, правило Марковникова)", "H2SO4 или H3PO4, t, p", [(1, a), (1, b)], [(1, S("alcohol", n, q))])
        if k == "KMnO4": ad("мягкое окисление (реакция Вагнера)", "H2O, 20 °C, раствор обесцвечивается", [(3, a), (2, b), (4, I("H2O"))], [(3, S("diol", n, p)), (2, I("MnO2")), (2, I("KOH"))])
        if b is None and p == 1 and n in (2, 3):
            ad("полимеризация", "t, p, катализатор", [("n ", a)], [(1, I("(-CH2-CH2-)n" if n == 2 else "(-CH2-CH(CH3)-)n"))])
    if c == "alkyne":
        if k == "H2":
            ad("гидрирование (неполное)", "Pd, катализатор Линдлара", [(1, a), (1, b)], [(1, S("alkene", n, p))])
            ad("гидрирование (полное)", "Ni, t", [(1, a), (2, b)], [(1, S("alkane", n))])
        if k == "HCl" and n == 2: ad("гидрохлорирование", "HgCl2, t", [(1, a), (1, b)], [(1, I("CH2=CH-Cl"))])
        if k == "H2O" and p == 1:
            ad("гидратация (реакция Кучерова)", "HgSO4, H2SO4, t", [(1, a), (1, b)], [(1, S("aldehyde", 2) if n == 2 else S("ketone", n, 2))])
        if b is None and n == 2: ad("тримеризация (реакция Зелинского)", "C активированный, 600 °C", [(3, a)], [(1, I("C6H6"))])
    if c == "alcohol":
        if k == "Na" and p == 1: ad("замещение (с активным металлом)", "", [(2, a), (2, b)], [(2, S("alkoxide", n)), (1, I("H2"))])
        if k in ("HCl", "HBr"): ad("замещение OH на галоген", "H2SO4 (конц.), t", [(1, a), (1, b)], [(1, S("hal" + X, n, p)), (1, I("H2O"))])
        if k == "H2SO4" and n > 1: ad("внутримолекулярная дегидратация (правило Зайцева)", "H2SO4 (конц.), t > 140 °C", [(1, a)], [(1, S("alkene", n, p)), (1, I("H2O"))])
        if k == "CuO":
            ad("окисление", "t", [(1, a), (1, b)], [(1, S("aldehyde", n) if p == 1 else S("ketone", n, p)), (1, I("Cu")), (1, I("H2O"))])
        if b and b[0] == "o" and b[1] == "acid" and p == 1:
            ad("этерификация", "H2SO4 (конц.), t", [(1, b), (1, a)], [(1, S("ester", b[2], n)), (1, I("H2O"))])
    if c == "aldehyde":
        if k == "H2": ad("присоединение (гидрирование)", "Ni, t", [(1, a), (1, b)], [(1, S("alcohol", n, 1))])
        if k == "Ag2O": ad("«серебряное зеркало»", "NH3 (р-р), t", [(1, a), (1, b)], [(1, S("acid", n)), (2, I("Ag"))])
        if k == "Cu(OH)2": ad("окисление свежим Cu(OH)2", "t", [(1, a), (2, b)], [(1, S("acid", n)), (1, I("Cu2O")), (2, I("H2O"))])
    if c == "ketone" and k == "H2": ad("присоединение (гидрирование)", "Ni, t", [(1, a), (1, b)], [(1, S("alcohol", n, p))])
    if c == "acid":
        if k == "NaOH": ad("нейтрализация", "", [(1, a), (1, b)], [(1, S("salt", n)), (1, I("H2O"))])
        if k == "Na2CO3": ad("обмен (выделение CO2)", "", [(2, a), (1, b)], [(2, S("salt", n)), (1, I("CO2")), (1, I("H2O"))])
        if k == "NaHCO3": ad("обмен (выделение CO2)", "", [(1, a), (1, b)], [(1, S("salt", n)), (1, I("CO2")), (1, I("H2O"))])
        if k == "Na": ad("замещение (с активным металлом)", "", [(2, a), (2, b)], [(2, S("salt", n)), (1, I("H2"))])
    if c.startswith("hal"):
        Y = c[3:]
        if k == "NaOH": ad("замещение (щелочной гидролиз)", "H2O, t", [(1, a), (1, b)], [(1, S("alcohol", n, p)), (1, I("Na" + Y))])
        if k == "KOH" and n > 1:
            ad("отщепление (дегидрогалогенирование)", "спиртовой раствор, t", [(1, a), (1, b)], [(1, S("alkene", n, p)), (1, I("K" + Y)), (1, I("H2O"))])
        if k == "Na" and p == 1 and 2 * n <= 6:
            ad("реакция Вюрца", "сухой эфир", [(2, a), (2, b)], [(1, S("alkane", 2 * n)), (2, I("Na" + Y))])
    return R


def pair(a, b):
    if a[0] == "i" and b[0] == "i":
        r = RXD.get(frozenset((a[1], b[1])))
        return [mk(*r)] if r else []
    out = []
    if a[0] == "o": out += org(a, b)
    if b[0] == "o": out += org(b, a)
    return out


def react(subs):
    res = []
    if len(subs) == 1 and subs[0][0] == "o": res += org(subs[0], None)
    for i in range(len(subs)):
        for j in range(i + 1, len(subs)):
            if subs[i] != subs[j]: res += pair(subs[i], subs[j])
    seen, out = set(), []
    for r in res:
        if r and r["eq"] not in seen:
            seen.add(r["eq"]); out.append(r)
    return out


def parse(raw):
    raw = raw.strip()
    if not raw: return None
    q = qn(raw)
    if q in NAME: return NAME[q]
    t = raw.translate(SUBT)
    if not re.search(r"[а-яё]", t): t = t.translate(LOOK)
    k = nk(t)
    return FORM.get(k) or FORMLOW.get(k.lower())


@app.route("/api/react", methods=["POST"])
def api():
    items = [str(x) for x in (request.get_json(silent=True) or {}).get("items", []) if str(x).strip()]
    subs, known, unknown = [], [], []
    for x in items:
        s = parse(x)
        if s is None: unknown.append(x); continue
        if s not in subs:
            subs.append(s); known.append({"f": INF[s][0], "n": INF[s][1]})
    return jsonify(known=known, unknown=unknown, reactions=react(subs) if subs else [])


@app.route("/")
def index():
    return Response(HTML, mimetype="text/html")


HTML = r'''<!doctype html>
<html lang="ru"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Химический конструктор реакций</title>
<style>
:root{--bg:#0a0c0b;--s1:#121614;--s2:#1a201d;--ln:#27322c;--tx:#e8f0eb;--mu:#8ea397;--g:#38c172}
*{box-sizing:border-box}
body{margin:0;background:var(--bg);color:var(--tx);font:16px/1.5 "Segoe UI",system-ui,sans-serif}
.w{max-width:860px;margin:0 auto;padding:22px 16px 70px}
h1{font-size:27px;margin:6px 0 4px;letter-spacing:-.3px}
.sub{color:var(--mu);margin:0 0 20px}
.card{background:var(--s1);border:1px solid var(--ln);border-radius:12px;padding:16px;margin-bottom:14px}
.card h2{font-size:15px;margin:0 0 10px;color:var(--mu);font-weight:600}
.row{display:flex;gap:8px;margin-bottom:8px}
input{flex:1;min-width:0;background:var(--s2);border:1px solid var(--ln);color:var(--tx);border-radius:8px;padding:11px 12px;font:17px Consolas,"JetBrains Mono",monospace;outline:none}
input:focus{border-color:var(--g)}
button{font:inherit;cursor:pointer;border-radius:8px;border:1px solid var(--ln);background:var(--s2);color:var(--tx);padding:8px 12px}
button:hover{border-color:var(--g)}
.pri{background:var(--g);border-color:var(--g);color:#06140c;font-weight:700}
.pri:hover{background:#4fd286}
.kb{display:flex;flex-wrap:wrap;gap:6px;margin-bottom:10px}
.kb button{font-family:Consolas,monospace;padding:7px 11px;min-width:42px}
.kb .rad{color:var(--g)}
.kt{font-size:13px;color:var(--mu);margin:8px 0 6px}
.acts{display:flex;gap:8px;flex-wrap:wrap;margin-top:6px}
.chips{display:flex;flex-wrap:wrap;gap:6px}
.chips button{font-size:14px;padding:6px 10px}
.f{font-family:Consolas,"JetBrains Mono",monospace}
.rec{padding:3px 0}
.bad{border-color:#6b3a2e;color:#f0c9bd}
.tag{display:inline-block;font-size:13px;color:var(--g);border:1px solid var(--g);border-radius:6px;padding:1px 8px;margin-bottom:8px}
.eq{font:19px/1.6 Consolas,"JetBrains Mono",monospace;overflow-x:auto;padding:4px 0;white-space:nowrap}
.cond{color:var(--mu);font-size:14px;margin-bottom:8px}
.pr{padding:8px 0;border-top:1px solid var(--ln);display:flex;justify-content:space-between;gap:10px;align-items:center}
.pr span.n{color:var(--mu)}
.use{font-size:13px;padding:4px 9px;white-space:nowrap}
</style></head><body><div class="w">
<h1>Химический конструктор</h1>
<p class="sub">Органика и неорганика, 10–11 класс. Вводи формулы (CH2=CH-CH3) или названия (пропен, уксусная кислота).</p>

<div class="card"><h2>Исходные вещества</h2>
<div id="rows"></div>
<div class="acts"><button id="add">+ вещество</button><button class="pri" id="go">Получить продукты</button></div></div>

<div class="card"><h2>Химическая клавиатура</h2><div id="kb"></div>
<div class="acts"><button id="bs">⌫ стереть</button><button id="clr">очистить поле</button></div></div>

<div class="card"><h2>Примеры по темам</h2><div class="chips" id="ex"></div></div>
<div id="out"></div></div>
<script>
const $=s=>document.querySelector(s);
const esc=s=>String(s).replace(/[&<>]/g,c=>({"&":"&amp;","<":"&lt;",">":"&gt;"}[c]));
const fm=s=>esc(s).replace(/([A-Za-z\)])(\d+)/g,"$1<sub>$2</sub>");
let cur=null;
function row(v){const d=document.createElement("div");d.className="row";
 d.innerHTML='<input placeholder="формула или название" autocomplete="off" spellcheck="false"><button title="убрать">×</button>';
 const i=d.querySelector("input");i.value=v||"";i.onfocus=()=>cur=i;i.onkeydown=e=>{if(e.key==="Enter")go()};
 d.querySelector("button").onclick=()=>{d.remove();if(cur===i)cur=null};$("#rows").appendChild(d);return i}
function ins(t){if(!cur||!document.contains(cur))cur=document.querySelector("#rows input")||row("");
 const s=cur.selectionStart,e=cur.selectionEnd;cur.value=cur.value.slice(0,s)+t+cur.value.slice(e);
 cur.focus();cur.setSelectionRange(s+t.length,s+t.length)}
const KB=[["Элементы",["C","H","O","N","S","Cl","Br","Na","K","Cu","Ag","Mn"],""],
 ["Связи, скобки, цифры",["-","=","≡","(",")",..."0123456789"],""],
 ["Радикалы и группы",["CH3","C2H5","C3H7","C6H5","OH","NH2","COOH","CHO","NO2"],"rad"]];
KB.forEach(([t,a,c])=>{$("#kb").insertAdjacentHTML("beforeend",'<div class="kt">'+t+'</div>');
 const d=document.createElement("div");d.className="kb";
 a.forEach(x=>{const b=document.createElement("button");b.textContent=x;b.className=c;
  b.onmousedown=e=>e.preventDefault();b.onclick=()=>ins(x);d.appendChild(b)});$("#kb").appendChild(d)});
$("#bs").onmousedown=e=>e.preventDefault();
$("#bs").onclick=()=>{if(!cur)return;const s=cur.selectionStart,e=cur.selectionEnd;
 if(s!==e)cur.setRangeText("",s,e,"end");else if(s>0)cur.setRangeText("",s-1,s,"end");cur.focus()};
$("#clr").onclick=()=>{if(cur){cur.value="";cur.focus()}};
$("#add").onclick=()=>row("").focus();
const EX={"Алканы":[["CH4","Cl2"],["C3H8","O2"],["этан"]],
"Алкены":[["CH2=CH-CH3","HBr"],["этилен","H2O"],["пропен","KMnO4"],["CH2=CH2"]],
"Алкины":[["ацетилен","H2O"],["CaC2","H2O"],["пропин","H2"],["C2H2"]],
"Спирты":[["этанол","Na"],["этанол","CuO"],["CH3-CH2-OH","CH3COOH"],["пропанол-2","H2SO4"]],
"Альдегиды и кислоты":[["этаналь","Ag2O"],["этаналь","H2"],["уксусная кислота","NaOH"],["CH3COOH","Na2CO3"]],
"Галогеналканы":[["бромэтан","NaOH"],["CH3-CH2-Br","KOH"],["хлорметан","Na"]],
"Неорганика":[["HCl","NaOH"],["Na","H2O"],["H2SO4","NaOH"]]};
Object.entries(EX).forEach(([t,l])=>{$("#ex").insertAdjacentHTML("beforeend",'<div class="kt" style="width:100%;margin:6px 0 0">'+t+'</div>');
 l.forEach(x=>{const b=document.createElement("button");b.textContent=x.join(" + ");
  b.onclick=()=>{$("#rows").innerHTML="";x.forEach(v=>row(v));go()};$("#ex").appendChild(b)})});
async function go(){
 const items=[...document.querySelectorAll("#rows input")].map(i=>i.value).filter(v=>v.trim());
 if(!items.length){$("#out").innerHTML='<div class="card bad">Введите хотя бы одно вещество.</div>';return}
 const r=await fetch("/api/react",{method:"POST",headers:{"Content-Type":"application/json"},body:JSON.stringify({items})});
 show(await r.json())}
$("#go").onclick=go;
function show(d){let h="";
 if(d.known.length)h+='<div class="card"><h2>Распознано</h2>'+d.known.map(k=>`<div class="rec"><span class="f">${fm(k.f)}</span> — ${esc(k.n)}</div>`).join("")+"</div>";
 if(d.unknown.length)h+=`<div class="card bad">Не распознано: <b>${d.unknown.map(esc).join(", ")}</b>.<br>Проверь написание. Поддерживаются неразветвлённые цепи до C6 и основные неорганические вещества.</div>`;
 if(!d.reactions.length&&d.known.length)h+='<div class="card">Для этого набора реакции из программы 10–11 класса не найдены. Попробуй реагенты: H2, Br2, HBr, H2O, Na, NaOH, CuO, O2, KMnO4, Ag2O.</div>';
 d.reactions.forEach(r=>{h+=`<div class="card"><span class="tag">${esc(r.type)}</span>
 <div class="eq">${fm(r.eq)}</div>${r.cond?`<div class="cond">Условия: ${fm(r.cond)}</div>`:""}
 ${r.prods.map(p=>`<div class="pr"><div><span class="f">${fm(p.f)}</span><br><span class="n">${esc(p.n)}</span></div><button class="use" data-f="${esc(p.f)}">использовать</button></div>`).join("")}</div>`});
 $("#out").innerHTML=h}
document.addEventListener("click",e=>{if(e.target.classList.contains("use")){row(e.target.dataset.f);window.scrollTo({top:0,behavior:"smooth"})}});
row("");row("");
</script></body></html>'''

if __name__ == "__main__":
    print("Открой в браузере: http://127.0.0.1:5050")
    app.run(host="127.0.0.1", port=5050, debug=False, use_reloader=False)
