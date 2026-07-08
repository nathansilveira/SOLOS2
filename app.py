# -*- coding: utf-8 -*-
"""
Análise de Estabilidade de Taludes — App interativo
Trabalho de Mecânica dos Solos II (UFSC) — Slope/W

Rodar localmente:  streamlit run app.py
"""
import numpy as np
import matplotlib.pyplot as plt
import streamlit as st

from slope_engine import (SURFACE, BASE_Y, WATER_INTERM, WATER_SAT,
                          ground_y, search_critical)

st.set_page_config(page_title='Estabilidade de Taludes', page_icon='⛰️',
                   layout='wide')

C_FELL = '#4285F4'
C_BISH = '#EA4335'
C_MP = '#FBBC04'


def fmt(v):
    return f'{v:.3f}'.replace('.', ',')


@st.cache_data(show_spinner=False)
def rodar_analise(gamma, c, phi, condition, q1, q2):
    return search_critical(gamma, c, phi, condition, q1, q2)


def desenhar_secao(res, condition, q1, q2):
    fig, ax = plt.subplots(figsize=(10, 5))

    # solo
    poly_x = np.concatenate([SURFACE[:, 0], [89, 0]])
    poly_y = np.concatenate([SURFACE[:, 1], [BASE_Y, BASE_Y]])
    ax.fill(poly_x, poly_y, color='#D2B48C', alpha=0.6, zorder=1)
    ax.plot(SURFACE[:, 0], SURFACE[:, 1], color='#5C4033', lw=2, zorder=3)

    # nível d'água
    if condition == 'intermediario':
        ax.plot(WATER_INTERM[:, 0], WATER_INTERM[:, 1], 'b--', lw=1.8,
                label='Nível freático', zorder=3)
    elif condition == 'saturado':
        ax.plot(WATER_SAT[:, 0], WATER_SAT[:, 1], 'b--', lw=1.8,
                label='Nível freático', zorder=3)

    # cargas
    def setas(x0, x1, y, cor, rotulo):
        for xa in np.linspace(x0, x1, 6):
            ax.annotate('', xy=(xa, y), xytext=(xa, y + 3),
                        arrowprops=dict(arrowstyle='->', color=cor, lw=1.6))
        ax.text((x0 + x1) / 2, y + 3.6, rotulo, ha='center', fontsize=10,
                color=cor, fontweight='bold')

    if q1 > 0:
        setas(3.5, 8.5, 32, '#CC0000', f'q₁ = {q1:.0f} kN/m')
    if q2 > 0:
        setas(25, 30, 24, '#CC0000', f'q₂ = {q2:.0f} kN/m')

    # superfície crítica
    if res and np.isfinite(res.get('fs_bishop', np.inf)):
        xc, yc, R = res['xc'], res['yc'], res['R']
        x1, x2 = res['span']
        xs = np.linspace(x1, x2, 200)
        ys = yc - np.sqrt(np.maximum(R**2 - (xs - xc)**2, 0))
        ax.plot(xs, ys, color='red', lw=2.5, zorder=4,
                label='Superfície crítica')
        ax.plot(xc, yc, 'r+', markersize=12, markeredgewidth=2)

    ax.set_xlim(-3, 92)
    ax.set_ylim(BASE_Y - 2, 48)
    ax.set_aspect('equal')
    ax.set_xlabel('x (m)')
    ax.set_ylabel('y (m)')
    ax.grid(alpha=0.3)
    ax.legend(loc='upper right', fontsize=9)
    fig.tight_layout()
    return fig


# ============================ INTERFACE ============================
st.title('⛰️ Análise de Estabilidade de Taludes')
st.caption('Mecânica dos Solos II — Talude escalonado (3 × 8 m, 1V:1,5H) · '
           'Método das fatias com busca automática da superfície crítica')

with st.sidebar:
    st.header('⚙️ Parâmetros')

    st.subheader('Carregamento')
    modo = st.radio('Definir cargas por:', ['Matrícula (N)', 'Valores manuais'])
    if modo == 'Matrícula (N)':
        N = st.slider('N — último dígito da matrícula', 0, 9, 0)
        q1 = 70 + 2 * N
        q2 = 40 + 3 * N
        st.info(f'q₁ = 70 + 2·{N} = **{q1} kN/m**\n\n'
                f'q₂ = 40 + 3·{N} = **{q2} kN/m**')
    else:
        q1 = st.slider('q₁ — carga no topo (kN/m)', 0, 300, 70, 5)
        q2 = st.slider('q₂ — carga no platô (kN/m)', 0, 300, 40, 5)

    aplicar_carga = st.toggle('Aplicar carregamento', value=True)

    st.subheader('Nível freático')
    cond_nome = st.radio('Condição:', ['Seco', 'Intermediário', 'Saturado'])
    condition = {'Seco': 'seco', 'Intermediário': 'intermediario',
                 'Saturado': 'saturado'}[cond_nome]

    st.subheader('Solo (Mohr-Coulomb)')
    c = st.number_input("c' — coesão (kPa)", 0.0, 100.0, 13.0, 1.0)
    phi = st.number_input("φ' — ângulo de atrito (°)", 5.0, 45.0, 25.0, 1.0)
    gamma = 18.0 if condition == 'seco' else 20.0
    st.caption(f'γ adotado = {gamma:.0f} kN/m³ '
               f'({"natural" if condition == "seco" else "saturado"})')

q1_eff = q1 if aplicar_carga else 0
q2_eff = q2 if aplicar_carga else 0

with st.spinner('Buscando a superfície crítica...'):
    res = rodar_analise(gamma, c, phi, condition, q1_eff, q2_eff)

# ---------------- Resultados ----------------
if not np.isfinite(res.get('fs_bishop', np.inf)):
    st.error('Nenhuma superfície válida encontrada com esses parâmetros.')
    st.stop()

fs_f, fs_b = res['fs_fell'], res['fs_bishop']

col1, col2, col3 = st.columns(3)
col1.metric('Fellenius', fmt(fs_f))
col2.metric('Bishop Simplificado', fmt(fs_b))
col3.metric('Morgenstern-Price (≈ Bishop)', fmt(fs_b),
            help='Para superfícies circulares, o M-P resulta em valores '
                 'praticamente idênticos ao Bishop (diferença < 1%).')

fs_min = min(fs_f, fs_b)
if fs_min < 1.0:
    st.error(f'🔴 **FS = {fmt(fs_min)} < 1,0 — talude em ruptura** '
             'nessa condição.')
elif fs_min < 1.2:
    st.warning(f'🟠 FS = {fmt(fs_min)} — abaixo de todos os níveis mínimos '
               'de segurança (1,2 a 1,5).')
elif fs_min < 1.5:
    st.warning(f'🟡 FS = {fmt(fs_min)} — estável, mas não atende o nível '
               'alto de segurança (FS ≥ 1,5).')
else:
    st.success(f'🟢 FS = {fmt(fs_min)} — atende ao nível alto de segurança '
               '(FS ≥ 1,5).')

st.pyplot(desenhar_secao(res, condition, q1_eff, q2_eff))
st.caption(f"Superfície crítica: centro ({res['xc']:.1f}, {res['yc']:.1f}) m, "
           f"raio {res['R']:.1f} m")

# ---------------- Comparação dos 5 casos ----------------
st.divider()
st.subheader('📊 Comparação — 5 casos do trabalho')
st.caption('Calcula os 5 casos com as cargas definidas na barra lateral.')

if st.button('Calcular os 5 casos'):
    casos = [
        ('Caso 1 — Seco, sem carga', 18, 'seco', 0, 0),
        ('Caso 2 — Seco, com carga', 18, 'seco', q1, q2),
        ('Caso 3 — Intermediário, com carga', 20, 'intermediario', q1, q2),
        ('Caso 4 — Saturado, com carga', 20, 'saturado', q1, q2),
        ('Caso 5 — Saturado, sem carga', 20, 'saturado', 0, 0),
    ]
    nomes, fs_fell_l, fs_bish_l = [], [], []
    barra = st.progress(0.0)
    for i, (nome, g, cond, a, b) in enumerate(casos):
        r = rodar_analise(g, c, phi, cond, a, b)
        nomes.append(nome.split('—')[1].strip())
        fs_fell_l.append(r['fs_fell'])
        fs_bish_l.append(r['fs_bishop'])
        barra.progress((i + 1) / 5)
    barra.empty()

    fig, ax = plt.subplots(figsize=(10, 5))
    x = np.arange(5)
    w = 0.35
    b1 = ax.bar(x - w / 2, fs_fell_l, w, color=C_FELL, label='Fellenius')
    b2 = ax.bar(x + w / 2, fs_bish_l, w, color=C_BISH,
                label='Bishop / M-P')
    for bars in (b1, b2):
        for bar in bars:
            h = bar.get_height()
            ax.text(bar.get_x() + bar.get_width() / 2, h + 0.02, fmt(h),
                    ha='center', va='bottom', fontsize=8.5,
                    fontweight='bold')
    ax.axhline(1.0, color='black', ls='--', lw=1.2)
    ax.text(4.55, 1.02, 'FS = 1,0', fontsize=9, style='italic')
    ax.set_xticks(x)
    ax.set_xticklabels([n.replace(', ', '\n') for n in nomes], fontsize=9)
    ax.set_ylabel('Fator de Segurança (FS)')
    ax.set_ylim(0, max(fs_bish_l) * 1.25)
    ax.legend(frameon=False)
    ax.grid(axis='y', alpha=0.3)
    fig.tight_layout()
    st.pyplot(fig)

    st.table({
        'Caso': [c_[0] for c_ in casos],
        'Fellenius': [fmt(v) for v in fs_fell_l],
        'Bishop / M-P': [fmt(v) for v in fs_bish_l],
    })

with st.expander('📋 Resultados de referência obtidos no Slope/W'):
    st.table({
        'Caso': ['1 — Seco, sem carga', '2 — Seco, com carga',
                 '3 — Intermediário, com carga', '4 — Saturado, com carga',
                 '5 — Saturado, sem carga'],
        'Fellenius': ['1,491', '1,401', '1,135', '0,908', '0,946'],
        'Bishop': ['1,564', '1,468', '1,246', '0,918', '0,945'],
        'Morgenstern-Price': ['1,562', '1,467', '1,247', '0,924', '0,950'],
    })

st.divider()
st.caption('⚠️ Cálculo pelo método das fatias com superfícies circulares, '
           'implementado em Python. Validado contra o Slope/W: diferença '
           'inferior a 1% nos casos sem carregamento; nos demais casos podem '
           'ocorrer pequenas diferenças pela posição da superfície crítica '
           'encontrada na busca automática. Ferramenta acadêmica — não '
           'substitui software de projeto.')
