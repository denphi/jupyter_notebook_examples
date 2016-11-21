APP = main_clpt
MFILES = 
matlabversion = 7.12
DELETE_FILES = 
BIN = ../bin
RP_DIR = $(shell . /etc/environ.sh; use -e -r rappture; echo $${RAPPTURE_PATH})

$(APP): $(APP).m $(MFILES)
	mcc -v -C -m -R -nojvm -I $(RP_DIR)/lib/matlab -I Library $(APP).m $(MFILES)
	patchrunscript.sh run_$(APP).sh run_$(APP).sh 1

all:
	. /etc/environ.sh; use -e -r matlab-$(matlabversion); make $(APP)

install: all
	install --mode 0755 -D run_$(APP).sh $(BIN)
	install --mode 0755 -D $(APP) $(BIN)
	if [ -e $(APP).ctf ] ; \
		then mv $(APP).ctf $(BIN); \
	fi
	cd $(BIN); \
	if [ -e $(APP).ctf ] ; \
	then \
		. /etc/environ.sh; \
		use -e -r matlab-$(matlabversion); \
		extractCTF $(APP).ctf; \
	fi
	install --mode 0644 -D readme.txt $(BIN)/$(APP).readme

clean:
	rm -f $(APP) $(APP).ctf \
         $(APP).prj run_$(APP).sh \
         $(APP)_main.c $(APP)_mcc_component_data.c \
         mccExcludedFiles.log readme.txt \
         $(DELETE_FILES)

distclean: clean
	rm -rf $(BIN)/$(APP) \
          $(BIN)/$(APP).ctf \
          $(BIN)/$(APP).readme \
          $(BIN)/$(APP)_mcr \
          $(BIN)/run_$(APP).sh

