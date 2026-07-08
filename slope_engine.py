"""
Motor de cálculo de estabilidade de taludes — método das fatias
Superfícies circulares com busca automática (Fellenius e Bishop Simplificado)

Geometria e parâmetros do trabalho de Mecânica dos Solos II (UFSC)
"""
import numpy as np

GAMMA_W = 9.81  # kN/m3

# ---------------- Geometria do talude ----------------
# Superfície do terreno (pontos A até I)
SURFACE = np.array([
    [0.0, 32.0], [12.0, 32.0], [24.0, 24.0], [31.0, 24.0],
    [43.0, 16.0], [50.0, 16.0], [62.0, 8.0], [69.0, 8.0], [89.0, 8.0]
])
BASE_Y = -4.0  # fundo do modelo

# Linhas piezométricas
WATER_INTERM = np.array([[0, 20], [12, 20], [24, 18], [31, 18], [43, 15], [89, 15]])
WATER_SAT = np.array([[0, 30], [12, 30], [24, 23], [31, 23], [43, 15],
                      [50, 15], [62, 7], [89, 7]])


def ground_y(x):
    return np.interp(x, SURFACE[:, 0], SURFACE[:, 1])


def water_y(x, condition):
    if condition == 'seco':
        return np.full_like(np.atleast_1d(x), -1e9, dtype=float)
    table = WATER_INTERM if condition == 'intermediario' else WATER_SAT
    return np.interp(x, table[:, 0], table[:, 1])


def surcharge(x, q1, q2):
    """Pressão vertical (kPa) na superfície na posição x.
    q1: faixa de 5 m centrada no topo (platô A-B, centro x=6)  -> 3,5 a 8,5
    q2: faixa de 5 m centrada no platô intermediário superior (C-D, centro x=27,5) -> 25 a 30
    """
    q = np.zeros_like(np.atleast_1d(x), dtype=float)
    q = np.where((x >= 3.5) & (x <= 8.5), q + q1, q)
    q = np.where((x >= 25.0) & (x <= 30.0), q + q2, q)
    return q


def circle_intersections(xc, yc, R):
    """Encontra o trecho em x onde o círculo está abaixo do terreno."""
    xs = np.linspace(max(xc - R, SURFACE[0, 0]), min(xc + R, SURFACE[-1, 0]), 2000)
    inside = R**2 - (xs - xc)**2
    valid = inside > 0
    if valid.sum() < 10:
        return None
    ybase = np.full_like(xs, np.nan)
    ybase[valid] = yc - np.sqrt(inside[valid])
    below = valid & (ybase < ground_y(xs))
    if below.sum() < 10:
        return None
    # maior trecho contíguo
    idx = np.where(below)[0]
    splits = np.where(np.diff(idx) > 1)[0]
    segments = np.split(idx, splits + 1)
    seg = max(segments, key=len)
    x1, x2 = xs[seg[0]], xs[seg[-1]]
    if x2 - x1 < 3.0:
        return None
    # o círculo precisa aflorar na superfície nas DUAS extremidades
    # (profundidade ~0 nos limites do trecho), senão a superfície é inválida
    d1 = ground_y(x1) - (yc - np.sqrt(max(R**2 - (x1 - xc)**2, 0.0)))
    d2 = ground_y(x2) - (yc - np.sqrt(max(R**2 - (x2 - xc)**2, 0.0)))
    if d1 > 0.6 or d2 > 0.6:
        return None
    return x1, x2


def analyze_circle(xc, yc, R, gamma, c, phi_deg, condition, q1, q2, n_slices=60):
    """Retorna (FS_fellenius, FS_bishop) para um círculo tentativa, ou None."""
    inter = circle_intersections(xc, yc, R)
    if inter is None:
        return None
    x1, x2 = inter
    phi = np.radians(phi_deg)
    tanphi = np.tan(phi)

    xs = np.linspace(x1, x2, n_slices + 1)
    xm = 0.5 * (xs[:-1] + xs[1:])
    b = np.diff(xs)

    inside = R**2 - (xm - xc)**2
    if np.any(inside <= 0):
        return None
    ybase = yc - np.sqrt(inside)
    if np.min(ybase) < BASE_Y:      # não pode passar do fundo do modelo
        return None
    ytop = ground_y(xm)
    h = ytop - ybase
    if np.any(h < 0) or np.max(h) < 1.0:
        return None

    # convenção: alfa positivo onde a base desce no sentido do movimento
    # (talude desce para a direita -> fatias do lado da crista têm xm < xc)
    sin_a = (xc - xm) / R
    cos_a = np.sqrt(inside) / R
    alpha = np.arcsin(np.clip(sin_a, -1, 1))
    l = b / cos_a

    W = gamma * h * b + surcharge(xm, q1, q2) * b   # peso + sobrecarga

    yw = water_y(xm, condition)
    # lâmina d'água acima do terreno (água empoçada) atua como sobrecarga
    W = W + GAMMA_W * np.maximum(0.0, yw - ytop) * b
    u = np.maximum(0.0, GAMMA_W * (yw - ybase))

    driving = np.sum(W * sin_a)

    # Empuxo hidrostático horizontal da água empoçada (lâmina acima do terreno)
    # sobre as faces que delimitam a massa deslizante:
    #  - na saída (jusante): empuxo aponta contra o movimento -> estabilizante
    #  - na entrada (montante): empuxo a favor do movimento -> instabilizante
    y_exit = float(ground_y(x2))
    yw_exit = float(np.atleast_1d(water_y(x2, condition))[0])
    if yw_exit > y_exit:
        hw = yw_exit - y_exit
        P = 0.5 * GAMMA_W * hw**2
        y_p = y_exit + hw / 3.0
        driving -= P * (yc - y_p) / R
    y_ent = float(ground_y(x1))
    yw_ent = float(np.atleast_1d(water_y(x1, condition))[0])
    if yw_ent > y_ent:
        hw = yw_ent - y_ent
        P = 0.5 * GAMMA_W * hw**2
        y_p = y_ent + hw / 3.0
        driving += P * (yc - y_p) / R

    if driving <= 1.0:
        return None

    # ----- Fellenius (método ordinário, variante do Slope/W: N' = (W - u.b).cos a) -----
    resist_f = np.sum(c * l + np.maximum(W - u * b, 0.0) * cos_a * tanphi)
    fs_fell = resist_f / driving

    # ----- Bishop simplificado (iterativo) -----
    fs = fs_fell
    for _ in range(60):
        m_a = cos_a * (1.0 + np.tan(alpha) * tanphi / fs)
        m_a = np.maximum(m_a, 0.2)
        resist_b = np.sum((c * b + np.maximum(W - u * b, 0.0) * tanphi) / m_a)
        fs_new = resist_b / driving
        if abs(fs_new - fs) < 1e-5:
            fs = fs_new
            break
        fs = fs_new
    return fs_fell, fs, (x1, x2)


def search_critical(gamma, c, phi_deg, condition, q1, q2,
                    xc_range=(0, 60), yc_range=(30, 90), refine=True):
    """Busca automática da superfície crítica (grade de centros + raios)."""
    best = {'fs_bishop': np.inf}
    for xc in np.arange(xc_range[0], xc_range[1] + 0.1, 3.0):
        for yc in np.arange(yc_range[0], yc_range[1] + 0.1, 3.0):
            d_surf = yc - ground_y(np.array([xc]))[0]
            r_min = max(d_surf + 2.0, 5.0)
            r_max = yc - BASE_Y - 0.3
            if r_max <= r_min:
                continue
            for R in np.linspace(r_min, r_max, 14):
                res = analyze_circle(xc, yc, R, gamma, c, phi_deg, condition, q1, q2)
                if res is None:
                    continue
                fs_f, fs_b, span = res
                if fs_b < best['fs_bishop']:
                    best = dict(fs_fell=fs_f, fs_bishop=fs_b,
                                xc=xc, yc=yc, R=R, span=span)
    if refine and np.isfinite(best['fs_bishop']):
        xc0, yc0, R0 = best['xc'], best['yc'], best['R']
        for xc in np.arange(xc0 - 3, xc0 + 3.01, 0.75):
            for yc in np.arange(yc0 - 3, yc0 + 3.01, 0.75):
                for R in np.arange(max(R0 - 4, 5), R0 + 4.01, 0.75):
                    res = analyze_circle(xc, yc, R, gamma, c, phi_deg, condition, q1, q2)
                    if res is None:
                        continue
                    fs_f, fs_b, span = res
                    if fs_b < best['fs_bishop']:
                        best = dict(fs_fell=fs_f, fs_bishop=fs_b,
                                    xc=xc, yc=yc, R=R, span=span)
    return best


if __name__ == '__main__':
    # Validação contra os resultados do Slope/W (N do aluno desconhecido -> teste sem carga)
    print('Caso 1 — seco, sem carga (Slope/W: Fell 1,491 | Bishop 1,564)')
    r = search_critical(gamma=18, c=13, phi_deg=25, condition='seco', q1=0, q2=0)
    print(f"  Fellenius = {r['fs_fell']:.3f} | Bishop = {r['fs_bishop']:.3f} | "
          f"centro=({r['xc']:.1f},{r['yc']:.1f}) R={r['R']:.1f}")

    print('Caso 5 — saturado, sem carga (Slope/W: Fell 0,946 | Bishop 0,945)')
    r = search_critical(gamma=20, c=13, phi_deg=25, condition='saturado', q1=0, q2=0)
    print(f"  Fellenius = {r['fs_fell']:.3f} | Bishop = {r['fs_bishop']:.3f} | "
          f"centro=({r['xc']:.1f},{r['yc']:.1f}) R={r['R']:.1f}")
