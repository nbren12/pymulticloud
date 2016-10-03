add_library (minpack STATIC dogleg.f dpmpar.f enorm.f fdjac1.f hybrd.f hybrd1.f qform.f qrfac.f r1mpyq.f r1updt.f)

add_library (minpackdl SHARED dogleg.f dpmpar.f enorm.f fdjac1.f hybrd.f hybrd1.f qform.f qrfac.f r1mpyq.f r1updt.f)


