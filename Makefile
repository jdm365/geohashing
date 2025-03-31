all: install clean

install:
	uv pip uninstall geohashing
	uv pip install -n .

clean:
	rm -rf build dist *.egg-info .cache
