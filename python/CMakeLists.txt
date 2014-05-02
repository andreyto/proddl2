include(PythonInstaller)

configure_file("config.py.in" ${PY_SOURCE_PACKAGE_DIR}/${PY_PACKAGE_NAME}/config.py @ONLY)

install_python_package(
	DEPS "${PY_SOURCE_PACKAGE_DIR}/${PY_PACKAGE_NAME}/__init__.py"
	)

message("WRAPPER=${WRAPPER}")
add_test(NAME pytest COMMAND ${WRAPPER}
	py.test --pyargs ${PY_PACKAGE_NAME})
