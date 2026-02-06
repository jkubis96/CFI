.PHONY: format lint all

format:
	isort cfi
	black cfi
	isort tests
	black tests


lint:
	pylint --exit-zero cfi/CellFunctionality.py
	pylint --exit-zero cfi/CellGraph.py



all: format lint
