.PHONY: format lint all

format:
	isort cfi_toolkit
	black cfi_toolkit
	isort tests
	black tests


lint:
	pylint --exit-zero cfi_toolkit/CellFunctionality.py
	pylint --exit-zero cfi_toolkit/CellGraph.py



all: format lint
