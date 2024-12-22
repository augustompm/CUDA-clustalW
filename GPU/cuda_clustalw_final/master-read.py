#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Script: generate_project_report.py
----------------------------------
Gera um relatório ultra-completo do projeto CUDA-clustalW:

- Mapeia diretórios, assinaturas de função
- Mostra referências a "sm_13", "sm_86", NVCCFLAGS, etc., exibindo snippet de código
- Localiza APIs obsoletas (cudaThreadSynchronize, etc.) e exibe snippet
- Exibe TODO, FIXME, etc. com snippet
- Separa warnings e errors do 'make' e agrupa
- Gera contagem estatística de cada tipo de referência

Uso:
    ./generate_project_report.py
    # ou python3 generate_project_report.py

Requer Python 3.
"""

import os
import re
import sys
import subprocess
from collections import defaultdict, Counter

# --------------------------------------------------------------------
# CONFIGURAÇÕES GERAIS
# --------------------------------------------------------------------

DEFAULT_DIR = os.path.expanduser("~/CUDA-clustalW/GPU/cuda_clustalw_final")
LOG_FILE = "project_report.txt"

# Número de linhas de contexto antes/depois de cada match
CONTEXT_LINES = 3

# Regexes compilados para assinaturas de função
FUNC_REGEX = re.compile(
    r"(?:__global__|__device__|__host__|\bstatic\b|\binline\b|\bextern\b|\b__forceinline__)?"
    r"\s+"
    r"[a-zA-Z_][a-zA-Z0-9_<>\*\s:]*"
    r"\s+"
    r"[a-zA-Z_][a-zA-Z0-9_]*"
    r"\s*\([^)]*\)"
)

# Regex (compilado) para compute capabilities/flags
CUDA_ARCH_REGEX = re.compile(
    r"(sm_\d+|compute_\d+|gencode\s*|arch=compute_\d+|arch=sm_\d+|NVCCFLAGS|CUDAFLAGS|--gpu-architecture|--gpu-code)"
)

# Listas de strings para buscas literais
CUDA_OBSOLETE_APIS = [
    r"cudaThreadSynchronize",
    r"cudaThreadExit",
    r"__syncthreadsCount",
    r"cutil\.h",
    r"helper_cuda\.h",
    r"texture\<",
    r"BindTexture",
]

PROBLEM_MARKERS = [
    r"TODO",
    r"FIXME",
    r"DEPRECATED",
    r"AW:",
    r"eight[ -_]years",
]

# --------------------------------------------------------------------
# FUNÇÕES AUXILIARES
# --------------------------------------------------------------------

def check_cuda_env_info():
    """Retorna string com a versão do nvcc e variáveis de ambiente relevantes."""
    info = []
    # nvcc version
    try:
        out = subprocess.check_output(["nvcc", "--version"], stderr=subprocess.STDOUT)
        info.append("nvcc version:\n" + out.decode("utf-8"))
    except Exception as e:
        info.append(f"Falha ao chamar 'nvcc --version': {e}")

    # Ambiente
    relevant_vars = ["CUDA_HOME", "CUDA_PATH", "PATH", "LD_LIBRARY_PATH"]
    for var in relevant_vars:
        val = os.environ.get(var, "Not set")
        info.append(f"{var}: {val}")
    return "\n".join(info)


def get_lines_with_context(filepath, line_num, context=CONTEXT_LINES):
    """
    Retorna as linhas [start_line, end_line] do arquivo, para mostrar contexto ao redor do 'line_num'.
    O 'context' é quantas linhas antes/depois exibir.
    """
    snippet = []
    try:
        with open(filepath, "r", encoding="utf-8", errors="ignore") as f:
            all_lines = f.readlines()

        # Indices em Python são 0-based, mas line_num possivelmente 1-based. Vamos assumir 1-based
        idx = line_num - 1
        start = max(0, idx - context)
        end = min(len(all_lines), idx + context + 1)
        for i in range(start, end):
            # Exemplo: "[ 123]   código..."
            snippet.append(f"[{i+1:4}] {all_lines[i].rstrip()}")
    except:
        # Se der erro, retornamos snippet vazio
        pass

    return snippet


def find_patterns_in_file(filepath, pattern_or_list):
    """
    Encontra todas as ocorrências de 'pattern_or_list' no arquivo. 
    Retorna lista de tuplas (line_num, matched_text).
      - line_num é inteiro (base 1)
      - matched_text é a string que casou
    Se 'pattern_or_list' for re.Pattern, usamos .finditer()
    Se for list[str], cada str é compilada e buscamos matches literais (ou regex).
    """
    results = []
    try:
        with open(filepath, "r", encoding="utf-8", errors="ignore") as f:
            content = f.read()

        # Precisamos também do offset para line number
        lines_offsets = []
        offset = 0
        for line_idx, line in enumerate(content.splitlines(True), start=1):
            lines_offsets.append((line_idx, offset))
            offset += len(line)

        # Função local para descobrir line number a partir de position:
        def find_line_num(pos):
            # Faz uma busca binária simples ou linear
            # mas como n pode ser grande, mas iremos simplificar
            ln = 1
            for (lidx, off) in reversed(lines_offsets):
                if pos >= off:
                    ln = lidx
                    break
            return ln

        if isinstance(pattern_or_list, re.Pattern):
            for match in pattern_or_list.finditer(content):
                matched_text = match.group(0)
                pos = match.start()
                line_num = find_line_num(pos)
                results.append((line_num, matched_text))
        else:
            # pattern_or_list é list de strings (regex strings).
            # Para cada pattern, compilamos e rodamos finditer
            for pat in pattern_or_list:
                pat_compiled = re.compile(pat)
                for match in pat_compiled.finditer(content):
                    matched_text = match.group(0)
                    pos = match.start()
                    line_num = find_line_num(pos)
                    results.append((line_num, matched_text))

    except Exception as e:
        print(f"[AVISO] Erro lendo {filepath}: {e}", file=sys.stderr)
    
    return results


def parse_functions_in_file(filepath):
    """Encontra assinaturas de função com FUNC_REGEX; retorna lista (line_num, func_signature)."""
    results = []
    try:
        with open(filepath, "r", encoding="utf-8", errors="ignore") as f:
            content = f.read()

        # Precisamos do line offset se quisermos exibir snippet. Mas vamos usar approach simples:
        # rodar finditer para achar o match e descobrir a line a partir do index
        lines_offsets = []
        offset = 0
        for line_idx, line in enumerate(content.splitlines(True), start=1):
            lines_offsets.append((line_idx, offset))
            offset += len(line)

        def find_line_num(pos):
            ln = 1
            for (lidx, off) in reversed(lines_offsets):
                if pos >= off:
                    ln = lidx
                    break
            return ln

        for match in FUNC_REGEX.finditer(content):
            sig = match.group(0)
            pos = match.start()
            ln = find_line_num(pos)
            results.append((ln, sig))
    except Exception as e:
        print(f"[AVISO] Erro lendo {filepath}: {e}", file=sys.stderr)

    return results


def walk_directory(root_dir, log_lines, summary_counts):
    """
    Percorre o diretório. Para cada arquivo, busca:
     - assinaturas de função
     - references a sm_13/sm_86 (CUDA_ARCH_REGEX)
     - references a APIs obsoletas (CUDA_OBSOLETE_APIS)
     - references a marcadores (PROBLEM_MARKERS)
     Exibe snippet ao redor e registra contadores em summary_counts.
    """
    code_exts = [".h", ".hpp", ".c", ".cpp", ".cu", ".cuh"]

    for dirpath, dirnames, filenames in os.walk(root_dir):
        rel_path = os.path.relpath(dirpath, root_dir)
        if rel_path == ".":
            rel_path = root_dir
        log_lines.append(f"\n[DIRECTORY] {rel_path}")

        for filename in filenames:
            filepath = os.path.join(dirpath, filename)
            log_lines.append(f"  |- {filename}")

            if not os.path.isfile(filepath):
                continue

            # 1) Se for extensão de código, checar assinaturas de função
            ext = os.path.splitext(filename)[1].lower()
            if ext in code_exts:
                funcs = parse_functions_in_file(filepath)
                if funcs:
                    log_lines.append(f"    -> Found {len(funcs)} function(s):")
                    for (ln, signature) in funcs:
                        log_lines.append(f"       [line {ln}] {signature}")
                        # Atualiza contagem
                        summary_counts["FUNCTIONS"] += 1

            # 2) Buscas de referência
            arch_refs = find_patterns_in_file(filepath, CUDA_ARCH_REGEX)
            obs_refs  = find_patterns_in_file(filepath, CUDA_OBSOLETE_APIS)
            prob_refs = find_patterns_in_file(filepath, PROBLEM_MARKERS)

            all_matches = []
            # arch_refs => label "CUDA Arch / Flags"
            for (ln, text) in arch_refs:
                all_matches.append(("CUDA Arch / Flags", ln, text))
                summary_counts["CUDA Arch / Flags"] += 1
            # obs_refs => label "API CUDA Obsoleta"
            for (ln, text) in obs_refs:
                all_matches.append(("API CUDA Obsoleta", ln, text))
                summary_counts["API CUDA Obsoleta"] += 1
            # prob_refs => label "Marcador de problema"
            for (ln, text) in prob_refs:
                all_matches.append(("Marcador de problema", ln, text))
                summary_counts["Marcador de problema"] += 1

            # Se existirem matches, printar snippet
            if all_matches:
                log_lines.append("    -> Referências importantes (com snippet):")
                for (lbl, ln, matched) in all_matches:
                    log_lines.append(f"       [{lbl}] line {ln}: '{matched}'")
                    snippet_lines = get_lines_with_context(filepath, ln, CONTEXT_LINES)
                    for snippet_line in snippet_lines:
                        log_lines.append(f"         {snippet_line}")


def run_make_and_capture(build_dir):
    """Roda make e separa warnings/erros para posterior parse e exibição."""
    try:
        proc = subprocess.Popen(
            ["make"], cwd=build_dir, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
        )
        out, err = proc.communicate()
        code = proc.returncode
    except FileNotFoundError:
        return -1, ("ERRO: 'make' não encontrado.\n", "")
    except Exception as e:
        return -2, (f"ERRO inesperado ao chamar make: {e}\n", "")
    
    return code, (out, err)


def parse_make_output(make_out, make_err):
    """
    Faz parse básico de warnings e errors em make_out/make_err.
    Retorna (all_warnings, all_errors) como listas de strings e
    (parsed_lines) com todo o output raw.
    """
    # Um parse simples: procurar linhas que contenham "error:" ou "warning:".
    # Exemplo: "cudaFullPairwiseAlign.cu(83): error: no instance..."
    # ou clang/gcc style: "file.cpp:123: error: something"
    warnings = []
    errors = []
    parsed_lines = []

    def classify_line(line):
        lower = line.lower()
        if "error:" in lower:
            errors.append(line)
        elif "warning:" in lower:
            warnings.append(line)

    # Processar stdout
    for line in make_out.splitlines():
        parsed_lines.append(line)
        classify_line(line)

    # Processar stderr
    for line in make_err.splitlines():
        parsed_lines.append(line)
        classify_line(line)

    return warnings, errors, parsed_lines


# --------------------------------------------------------------------
# MAIN
# --------------------------------------------------------------------

def main():
    log_lines = []
    # Cabeçalho
    header = (
        "============================================================\n"
        " Relatório Ultra-Completo do Projeto CUDA-clustalW\n"
        "============================================================\n\n"
        "Inclui:\n"
        "1) Árvore de diretórios, assinaturas de função\n"
        "2) Snippets para menções a sm_13, sm_86, APIs obsoletas, TODO, etc.\n"
        "3) Estatísticas ao final\n"
        "4) Execução do 'make', com parse de warnings e errors\n"
        "5) Enfim, informações detalhadas para migração/debug.\n\n"
    )
    log_lines.append(header)

    # Checa ambiente
    log_lines.append("=== Informações de ambiente ===\n")
    log_lines.append(check_cuda_env_info())
    log_lines.append("\n")

    # Sumário (contadores) para final do relatório
    summary_counts = Counter()  # ex.: summary_counts["CUDA Arch / Flags"] += 1

    # Caminho do projeto
    target_dir = DEFAULT_DIR
    if not os.path.isdir(target_dir):
        log_lines.append(f"ERRO: Diretório '{target_dir}' não encontrado.\n")
    else:
        # Mapeia
        log_lines.append(f"=== Mapeando o diretório: {target_dir} ===\n")
        walk_directory(target_dir, log_lines, summary_counts)

        # Rodar make
        log_lines.append("\n=== Execução do 'make' ===\n")
        code, (mk_out, mk_err) = run_make_and_capture(target_dir)
        log_lines.append(f"Código de retorno do 'make': {code}\n")

        # Parse warnings e errors
        wlist, elist, raw_lines = parse_make_output(mk_out, mk_err)

        log_lines.append(">> RAW MAKE OUTPUT <<\n")
        if not mk_out.strip() and not mk_err.strip():
            log_lines.append("(nenhuma saída do make)\n")
        else:
            for l in raw_lines:
                log_lines.append(l)

        # Exibe warnings e errors em seções separadas
        log_lines.append("\n=== Warnings Coletados ===\n")
        if wlist:
            summary_counts["MAKE_WARNINGS"] = len(wlist)
            for w in wlist:
                log_lines.append("- " + w)
        else:
            log_lines.append("(nenhum warning encontrado)")

        log_lines.append("\n=== Errors Coletados ===\n")
        if elist:
            summary_counts["MAKE_ERRORS"] = len(elist)
            for e in elist:
                log_lines.append("* " + e)
        else:
            log_lines.append("(nenhum erro encontrado)")

    # Sumário final
    log_lines.append("\n\n=== SUMÁRIO DE CONTAGEM DE ITENS DETECTADOS ===\n")
    for key, val in summary_counts.most_common():
        log_lines.append(f"- {key}: {val}")

    # Salvar
    with open(LOG_FILE, "w", encoding="utf-8") as f:
        f.write("\n".join(log_lines))

    print(f"Relatório gerado em: {LOG_FILE}")


if __name__ == "__main__":
    main()
