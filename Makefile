.DEFAULT_GOAL := test

#
# Virtual environment setup
# You may want to replace python3 with the path to your python3 executable,
# e.g. the output of `pyenv which python3`when using pyenv
#
create-venv:
	python3 -m venv '.venv'
	echo "Don't forget to activate with 'source .venv/bin/activate'"

#
# Install requirements
#
install:
	python3 -m pip install --upgrade pip
	python3 -m pip install flit poetry
	python3 -m flit install --deps production

install-dev: install
	python3 -m pip install '.[dev]'
	python3 -m certifi
	pre-commit install --hook-type pre-commit --hook-type pre-push
	echo "Please also install pandoc to create the documentation."

#
# Checks & package upload
#
check: install-dev
	isort pytximport
	black --check pytximport
	flake8 pytximport
	bandit -ll --recursive pytximport

upload: check
	flit publish

#
# Testing
#
unittest:
	coverage run -m pytest --maxfail=10

coverage-report: unittest
	coverage report
	coverage html

test: check unittest
