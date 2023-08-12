(asdf:defsystem convex-covering
  :version "1.0.0"
  :license "zlib"
  :author "Yukari Hafner <shinmera@tymoon.eu>"
  :maintainer "Yukari Hafner <shinmera@tymoon.eu>"
  :description "An algorithm for the computation of a convex hull covering set of a mesh"
  :homepage "https://shirakumo.github.io/convex-covering/"
  :bug-tracker "https://github.com/shirakumo/convex-covering/issues"
  :source-control (:git "https://github.com/shirakumo/convex-covering.git")
  :serial T
  :components ((:file "package")
               (:file "decomposition")
               (:file "documentation"))
  :depends-on (:manifolds
               :3d-spaces
               :quickhull
               :priority-queue
               :documentation-utils)
  :in-order-to ((asdf:test-op (asdf:test-op :convex-covering-test))))

