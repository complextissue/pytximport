.DEFAULT_GOAL := test

#
# Install requirements
#
install:
	uv sync

install-dev:
	uv sync --all-extras
	pre-commit install --hook-type pre-commit --hook-type pre-push
	echo "Please also install pandoc to create the documentation."

#
# Checks
#
check: install-dev
	uv run ruff format --check pytximport
	uv run ruff check pytximport
	uv run mypy -p pytximport
	uv run bandit -ll --recursive pytximport

#
# Testing
#
unittest:
	uv run coverage run -m pytest --maxfail=10

coverage-report: unittest
	uv run coverage report
	uv run coverage html

test: check unittest
