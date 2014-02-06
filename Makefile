all: RsplineCov.so RmixCov.so

RsplineCov.so: src/splineCov.cpp
	R CMD SHLIB -o RsplineCov.so src/splineCov.cpp

RmixCov.so: src/mixCov.cpp
	R CMD SHLIB -o RmixCov.so src/mixCov.cpp

clean:
	rm -f src/*.o
	rm -f RsplineCov.so
	rm -f RmixCov.so
