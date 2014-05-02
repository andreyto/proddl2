#!/bin/sh
g++ -O2 -I../../../include -DDBG_ENABLED -DATLOG_LEVEL=5 -o DETest DETest.cpp DESolver.cpp ../External/Dbg/dbg.cpp ../Common/logger.cpp

