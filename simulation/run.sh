#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SRC_DIR="$SCRIPT_DIR/src/main/java"
BIN_DIR="$SCRIPT_DIR/bin"
MAIN_CLASS="ar.edu.itba.sds.tp3.simulation.SimulationMain"
JAVAC_CMD="javac"
JAVA_CMD="java"

if ! command -v javac >/dev/null 2>&1; then
	if [[ -d "$HOME/.jdks" ]]; then
		FALLBACK_JAVAC="$(find "$HOME/.jdks" -maxdepth 4 -type f -name "javac.exe" 2>/dev/null | head -n 1)"
		if [[ -n "$FALLBACK_JAVAC" ]]; then
			JAVAC_CMD="$FALLBACK_JAVAC"
			JAVA_CMD="${FALLBACK_JAVAC%/javac.exe}/java.exe"
		fi
	fi

	if [[ "$JAVAC_CMD" == "javac" ]]; then
		echo "Error: javac was not found in PATH. Install a JDK (17+) and retry." >&2
		exit 1
	fi
fi

mkdir -p "$BIN_DIR"

mapfile -t JAVA_FILES < <(find "$SRC_DIR" -name "*.java" | sort)

if [[ ${#JAVA_FILES[@]} -eq 0 ]]; then
	echo "Error: no Java sources were found under $SRC_DIR" >&2
	exit 1
fi

JAVAC_OUTPUT_DIR="$BIN_DIR"
JAVA_CLASSPATH_DIR="$BIN_DIR"

if [[ "$JAVAC_CMD" == *.exe ]]; then
	JAVAC_OUTPUT_DIR="$(cygpath -w "$BIN_DIR")"
	JAVA_CLASSPATH_DIR="$(cygpath -w "$BIN_DIR")"
fi

"$JAVAC_CMD" -d "$JAVAC_OUTPUT_DIR" "${JAVA_FILES[@]}"
"$JAVA_CMD" -cp "$JAVA_CLASSPATH_DIR" "$MAIN_CLASS" "$@"