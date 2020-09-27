## change the following setting for different machines 
GEOTIFFINC=-I/opt/Toolkit-x86_64/libgeotiff/include
TIFFINC=-I/opt/Toolkit-x86_64/libtiff-4.0.6/include
GEOTIFFLIB=-L/opt/Toolkit-x86_64/libgeotiff/lib/ -lgeotiff  
   TIFFLIB=-L/opt/Toolkit-x86_64/libtiff-4.0.6/lib -ltiff 
GEOTIFFLIB=-L/opt/Toolkit-x86_64/libgeotiff/lib/ -lgeotiff  
   TIFFLIB=-L/opt/Toolkit-x86_64/libtiff-4.0.6/lib -ltiff 
   GCTPLIB=-L/opt/Toolkit-x86_64/hdfeos2-19.1.00/lib -lGctp 


CC=gcc

CPPFLAGS=-DMAX_NC_NAME=H4_MAX_NC_NAME -DMAX_VAR_DIMS=H4_MAX_VAR_DIMS -Wno-write-strings ${TIFFINC} ${GEOTIFFINC}

TARGET = Landsat.BRDF.normalize.v1.0
  OBJS = $(TARGET).o \
		 model.o \
		 ard.tif.io.o \
		 local_solar.o 

$(TARGET): $(OBJS)
	$(CC) ${CFLAGS} -o $(TARGET) $(OBJS) ${GEOTIFFLIB} ${TIFFLIB} ${GCTPLIB} -lm 

$(TARGET).o: $(TARGET).c
	$(CC) ${CFLAGS} -c $(TARGET).c ${CPPFLAGS}

model.o: model.c
	$(CC) $(CFLAGS) -c model.c ${CPPFLAGS}

ard.tif.io.o: ard.tif.io.c
	$(CC) $(CFLAGS) -c ard.tif.io.c ${CPPFLAGS}

local_solar.o: local_solar.c
	$(CC) $(CFLAGS) -c local_solar.c ${CPPFLAGS}

clean:
	rm -f *.o

delete:
	rm -f $(TARGET)

#$(TARGET): $(OBJS)
#	$(CC) ${CFLAGS} -o $(TARGET) $(OBJS) $(OPT) ${CPPFLAGS} -lm

#CPPFLAGS=-DMAX_NC_NAME=H4_MAX_NC_NAME -DMAX_VAR_DIMS=H4_MAX_VAR_DIMS -Wno-write-strings ${TIFFINC} ${GEOTIFFINC} ${GCTPINC} ${HDFINC} ${HDFEOSINC} ${WELDINC}
#HDFINC=-I/opt/Toolkit-x86_64/hdf-4.2.11/include
#HDFEOSINC=-I/opt/Toolkit-x86_64/hdfeos2-19.1.00/include
#HDFINC=-I/opt/Toolkit-x86_64/hdf-4.2.11/include
#GCTPINC=-I/opt/Toolkit-x86_64/hdfeos5-1.15/include
