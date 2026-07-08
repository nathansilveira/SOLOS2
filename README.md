# Análise de Estabilidade de Taludes — App Interativo

Site interativo do trabalho de Mecânica dos Solos II (Slope/W).
A pessoa escolhe a condição do nível freático e ajusta as cargas (q1 e q2),
e o app recalcula o fator de segurança pelo método das fatias
(Fellenius e Bishop Simplificado) com busca automática da superfície crítica.

## Rodar no seu computador

```
pip install -r requirements.txt
streamlit run app.py
```

O site abre em http://localhost:8501

## Publicar de graça (Streamlit Community Cloud)

1. Crie uma conta no GitHub (se não tiver) e crie um repositório novo;
2. Envie os 3 arquivos: `app.py`, `slope_engine.py` e `requirements.txt`;
3. Acesse https://share.streamlit.io e faça login com o GitHub;
4. Clique em "Create app", escolha o repositório e o arquivo `app.py`;
5. Pronto — o site fica público em um link do tipo
   `https://seuapp.streamlit.app` pra mandar pra qualquer pessoa.

## Arquivos

- `slope_engine.py` — motor de cálculo (geometria, poropressão, Fellenius e
  Bishop, busca da superfície crítica). Validado contra o Slope/W:
  diferença < 1% nos casos sem carga.
- `app.py` — interface do site (controles, seção transversal com a
  superfície crítica, comparação dos 5 casos).
