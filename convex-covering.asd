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

               (:file "parallelization")

               (:file "debug")

               (:file "geometry")
               (:file "context")
               (:file "structures")
               (:file "flat-mesh")
               (:file "merging")
               (:file "decomposition")
               (:file "documentation"))
  :depends-on ("manifolds"
               "3d-spaces"
               "quickhull"
               "damn-fast-priority-queue"
               (:feature :convex-covering-with-lparallel "lparallel")
               "machine-state"
               "documentation-utils")
  :in-order-to ((asdf:test-op (asdf:test-op "convex-covering/test"))))

(defsystem "convex-covering/test"
  :depends-on ("parachute"

               "convex-covering"
               "convex-covering/visualization")
  :components ((:file "test"))
  :perform (test-op (operation component)
             (uiop:symbol-call "ORG.SHIRAKUMO.FRAF.CONVEX-COVERING.TEST"
                               '#:decompose-test-files)))

(defsystem "convex-covering/visualization"
  :depends-on ("alexandria"
               "cl-wavefront"
               "cl-dot"
               "inferior-shell"

               "convex-covering")
  :serial T
  :components ((:file "geometry-debug")
               (:file "export-hulls")
               (:file "visualization")))
