(defpackage #:org.shirakumo.fraf.convex-covering
  (:use #:cl #:org.shirakumo.fraf.math)
  (:local-nicknames
   (#:manifolds #:org.shirakumo.fraf.manifolds))
  (:export
   #:convex-hull
   #:vertices
   #:faces
   #:decompose))
