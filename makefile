objects=bige_estimation_module.o bige_objectives_module.o bige_input_module.o bige_main.o bige_output_module.o bige_operator_module.o

headers=bige_estimation_module.h bige_objectives_module.h bige_output_module.h bige_input_module.h bige_operator_module.h bige_define.h bige_struct.h

bige: $(objects)
	gcc -o bige $(objects) -ldl *.so -lm

bige_estimation_module.o: bige_estimation_module.h bige_external.h bige_define.h bige_struct.h
bige_objectives_module.o: bige_objectives_module.h bige_external.h bige_struct.h
bige_input_module.o: bige_input_module.h bige_external.h 
bige_output_module.o: bige_output_module.h bige_struct.h bige_external.h
bige_main.o: $(headers)
bige_operator_module.o: bige_estimation_module.h bige_objectives_module.h bige_operator_module.h bige_define.h bige_struct.h

.PHONY: clean
clean:
	rm bige $(objects)
