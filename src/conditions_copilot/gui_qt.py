from __future__ import annotations
import os, json, subprocess, sys
from typing import Optional, Dict, Any

# Try PyQt6 first, fallback to PySide6
try:
    from PyQt6.QtWidgets import (
        QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
        QPlainTextEdit, QLineEdit, QPushButton, QFileDialog, QMessageBox,
        QInputDialog, QLabel
    )
    from PyQt6.QtCore import Qt
    QT_LIB = "PyQt6"
except Exception:
    from PySide6.QtWidgets import (
        QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
        QPlainTextEdit, QLineEdit, QPushButton, QFileDialog, QMessageBox,
        QInputDialog, QLabel
    )
    from PySide6.QtCore import Qt
    QT_LIB = "PySide6"

from .controller import build_discovery_payload, execute_actions
from .schemas import ProposalsRequest
from .tools import dicts as dicts_mod, validator
from .llm_prompting.client import call_openai


def _read_text(path: str) -> str:
    try:
        with open(path, "r", encoding="utf-8") as f:
            return f.read()
    except Exception as e:
        return f"<error reading {path}: {e}>"


def _json_extract(text: str) -> str:
    """Extract outermost JSON object if text contains extra logs."""
    first = text.find("{")
    last = text.rfind("}")
    return text[first:last + 1] if first >= 0 and last >= 0 else text


def call_llm_gui(payload: Dict[str, Any], system_path: str, llm_cmd: Optional[str], parent: QWidget) -> Dict[str, Any]:
    system_txt = _read_text(system_path)
    # OpenAI-compatible path if no shell cmd and API key present (OpenAI or DashScope)
    if not llm_cmd:
        model = os.environ.get("OPENAI_MODEL") or os.environ.get("DASHSCOPE_MODEL") or "deepseek-v3.1"
        base_url = (
            os.environ.get("OPENAI_BASE_URL")
            or os.environ.get("DASHSCOPE_BASE_URL")
            or "https://dashscope.aliyuncs.com/compatible-mode/v1"
        )
        mb = QMessageBox(parent)
        mb.setIcon(QMessageBox.Icon.Question)
        mb.setWindowTitle("Send to API?")
        mb.setText(f"Send payload to model '{model}'?")
        mb.setInformativeText("Set OPENAI_MODEL/OPENAI_BASE_URL or DASHSCOPE_MODEL/DASHSCOPE_BASE_URL via environment.")
        mb.setStandardButtons(QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No)
        if mb.exec() == QMessageBox.StandardButton.Yes:
            try:
                return call_openai(payload, system_txt, model=model, base_url=base_url, json_mode=True)
            except Exception as e:
                raise RuntimeError(str(e))
        # fall through to paste mode if declined

    if not llm_cmd:
        # Interactive paste dialog
        prompt_preview = json.dumps(payload, indent=2)
        msg = "Paste LLM JSON reply for the payload shown in the transcript."
        # Use a multi-line dialog to paste JSON
        text, ok = QInputDialog.getMultiLineText(parent, "LLM Reply (JSON)", msg, "")
        if not ok or not text.strip():
            raise RuntimeError("No JSON pasted from LLM.")
        try:
            return json.loads(_json_extract(text))
        except Exception as e:
            raise RuntimeError(f"Failed to parse LLM JSON: {e}")

    # Non-interactive: external command
    prompt = f"SYSTEM:\n{system_txt}\n\nUSER:\n{json.dumps(payload)}"
    proc = subprocess.run(llm_cmd, input=prompt.encode("utf-8"), shell=True,
                          stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if proc.returncode != 0:
        raise RuntimeError(f"LLM command failed: {proc.stderr.decode('utf-8', 'ignore')}")
    out = proc.stdout.decode("utf-8", "ignore").strip()
    try:
        return json.loads(_json_extract(out))
    except Exception as e:
        raise RuntimeError(f"Failed to parse LLM JSON: {e}\nRaw: {out[:4000]}")


class ChatWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle(f"Conditions Copilot (GUI, {QT_LIB})")

        # Defaults
        base_dir = os.path.dirname(__file__)
        self.dataset_csv = os.path.join(base_dir, "..", "..", "examples", "datasets", "demo_neighbors.csv")
        self.dict_dir = os.path.join(base_dir, "..", "..", "data", "dicts")
        self.system_path = os.path.join(base_dir, "llm_prompting", "system.txt")
        self.llm_cmd = os.environ.get("LLM_CMD", "")

        central = QWidget(self)
        self.setCentralWidget(central)
        layout = QVBoxLayout(central)

        # Transcript
        self.transcript = QPlainTextEdit(self)
        self.transcript.setReadOnly(True)
        self.transcript.setPlaceholderText("System responses and tool results will appear here…")
        layout.addWidget(self.transcript, stretch=1)

        # Status row (paths/LLM)
        status_row = QHBoxLayout()
        self.lbl_llm = QLabel(self)
        self._refresh_status()
        btn_dataset = QPushButton("Dataset…", self)
        btn_dataset.clicked.connect(self.choose_dataset)
        btn_dicts = QPushButton("Dictionaries…", self)
        btn_dicts.clicked.connect(self.choose_dictdir)
        btn_system = QPushButton("System…", self)
        btn_system.clicked.connect(self.choose_system)
        btn_llm = QPushButton("LLM Cmd…", self)
        btn_llm.clicked.connect(self.set_llm_cmd)
        for w in (self.lbl_llm, btn_dataset, btn_dicts, btn_system, btn_llm):
            status_row.addWidget(w)
        status_row.addStretch(1)
        layout.addLayout(status_row)

        # Input row
        input_row = QHBoxLayout()
        self.input = QLineEdit(self)
        self.input.setPlaceholderText("Enter reaction SMILES, then press Send…")
        self.btn_send = QPushButton("Send", self)
        self.btn_send.clicked.connect(self.on_send)
        input_row.addWidget(self.input, stretch=1)
        input_row.addWidget(self.btn_send)
        layout.addLayout(input_row)

        self.resize(900, 700)

    def _refresh_status(self):
        mode = "external" if self.llm_cmd else "paste"
        self.lbl_llm.setText(f"LLM: {mode}")

    def _append(self, role: str, content: str):
        self.transcript.appendPlainText(f"[{role}]\n{content}\n")
        self.transcript.verticalScrollBar().setValue(self.transcript.verticalScrollBar().maximum())

    def choose_dataset(self):
        path, _ = QFileDialog.getOpenFileName(self, "Select dataset CSV", os.path.dirname(self.dataset_csv), "CSV Files (*.csv);;All Files (*)")
        if path:
            self.dataset_csv = path
            self._append("INFO", f"Dataset set to: {path}")

    def choose_dictdir(self):
        path = QFileDialog.getExistingDirectory(self, "Select dictionaries directory", self.dict_dir)
        if path:
            self.dict_dir = path
            self._append("INFO", f"Dictionaries dir set to: {path}")

    def choose_system(self):
        path, _ = QFileDialog.getOpenFileName(self, "Select system.txt", os.path.dirname(self.system_path), "Text Files (*.txt);;All Files (*)")
        if path:
            self.system_path = path
            self._append("INFO", f"System prompt set to: {path}")

    def set_llm_cmd(self):
        text, ok = QInputDialog.getText(self, "LLM Command", "Shell command (empty = paste mode):", text=self.llm_cmd)
        if ok:
            self.llm_cmd = text.strip()
            self._refresh_status()
            mode = "external" if self.llm_cmd else "paste"
            self._append("INFO", f"LLM mode: {mode}")

    def on_send(self):
        rxn = self.input.text().strip()
        if not rxn:
            return
        self.input.clear()
        try:
            self._append("USER", f"RXN SMILES: {rxn}")

            # Round 1: discovery
            disc = build_discovery_payload(rxn)
            self._append("SYSTEM", f"Discovery payload:\n{json.dumps(disc, indent=2)}")
            reply1 = call_llm_gui(disc, self.system_path, self.llm_cmd, self)
            self._append("LLM", json.dumps(reply1, indent=2))

            actions = reply1.get("actions", [])
            tool_results = execute_actions(actions, rxn, self.dataset_csv, self.dict_dir)
            self._append("TOOLS", json.dumps(tool_results, indent=2))

            # Gather for proposals
            feats: Dict[str, Any] = {}
            retr = {"neighbors": [], "coverage": {"n_total": 0, "n_close": 0, "chemotype_overlap": "low", "median_yield_close": None}, "negatives": []}
            dicts = dicts_mod.load_dicts(self.dict_dir)
            for _, res in tool_results.items():
                if res.get("ok") and "features" in res:
                    feats = res["features"]
                if res.get("ok") and "neighbors" in res:
                    retr["neighbors"] = res["neighbors"]
                    retr["coverage"] = res["coverage"]
                if res.get("ok") and "dictionaries" in res:
                    dicts = res["dictionaries"]

            req = {
                "reaction": {
                    "smiles": rxn,
                    "type_hint": (reply1.get("classification", {}) or {}).get("name"),
                    "class_confidence": (reply1.get("classification", {}) or {}).get("confidence"),
                },
                "features": feats,
                "retrieval": retr,
                "dictionaries": dicts,
                "constraints": build_discovery_payload(rxn).get("constraints"),
                "ask": ("Propose top 3 ConditionCore and full conditions; include 2×3 screening panel if weak support. "
                         "Use only dictionary items. Return JSON only."),
            }
            req = ProposalsRequest(**req).model_dump(by_alias=True)
            self._append("SYSTEM", f"Proposals request:\n{json.dumps(req, indent=2)}")

            reply2 = call_llm_gui(req, self.system_path, self.llm_cmd, self)
            self._append("LLM", json.dumps(reply2, indent=2))

            report = validator.validate_proposals(reply2, dict_dir=self.dict_dir)
            self._append("VALIDATION", json.dumps(report, indent=2))

        except Exception as e:
            QMessageBox.critical(self, "Error", str(e))


def main():
    app = QApplication(sys.argv)
    w = ChatWindow()
    w.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
