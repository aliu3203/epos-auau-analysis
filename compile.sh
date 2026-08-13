cd $EPO
rm -Rf $LIBDIR && cmake -B$LIBDIR && make -C$LIBDIR -j8
