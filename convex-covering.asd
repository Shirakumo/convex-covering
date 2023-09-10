(asdf:defsystem "convex-covering"
  :version "1.0.0"
  :license "zlib"
  :author ("Yukari Hafner <shinmera@tymoon.eu>"
           "Jan Moringen <jmoringe@techfak.uni-bielefeld.de>")
  :maintainer "Yukari Hafner <shinmera@tymoon.eu>"
  :description "An algorithm for the computation of a convex hull covering set of a mesh"
  :homepage "https://shirakumo.github.io/convex-covering/"
  :bug-tracker "https://github.com/shirakumo/convex-covering/issues"
  :source-control (:git "https://github.com/shirakumo/convex-covering.git")
  :serial T
  :components ((:file "package")
               (:file "geometry")
               (:file "merging")
               (:file "decomposition")
               (:file "documentation"))
  :depends-on ("manifolds"
               "3d-spaces"
               "quickhull"
               "documentation-utils")
  :in-order-to ((asdf:test-op (asdf:test-op "convex-covering/test"))))

(defsystem "convex-covering/test"
  :depends-on ("convex-covering"
               "parachute"
               "cl-wavefront")
  :components ((:file "test"))
  :perform (test-op (operation component)
                    (uiop:symbol-call "ORG.SHIRAKUMO.FRAF.CONVEX-COVERING.TEST"
                                      '#:decompose-test-files)))

(defsystem "convex-covering/debug"
  :depends-on ("alexandria"
               "cl-dot"
               "inferior-shell"

               "convex-covering"
               "convex-covering/test")
  :components ((:file "debug")))
