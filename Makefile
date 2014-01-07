all: RsplineCov.so

RsplineCov.so: src/splineCov.cpp
	R CMD SHLIB -o RsplineCov.so src/splineCov.cpp

clean:
	rm -f src/*.o
	rm -f RsplineCov.so
