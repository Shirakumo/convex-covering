(defpackage #:org.shirakumo.fraf.convex-covering.test
  (:use #:cl #:parachute)
  (:local-nicknames
   (#:m #:org.shirakumo.fraf.math)
   (#:convex-covering #:org.shirakumo.fraf.convex-covering)
   (#:wavefront #:org.shirakumo.fraf.wavefront))
  (:shadow
   #:test)
  (:export
   #:decompose-file
   #:decompose-test-files))

(in-package #:org.shirakumo.fraf.convex-covering.test)

(defvar *here* #.(make-pathname :name NIL :type NIL :defaults (or *compile-file-pathname* *load-pathname* (error "COMPILE-FILE or LOAD this."))))

(define-test convex-covering)

(defun obj-file (file)
  (merge-pathnames (merge-pathnames file (make-pathname :type "obj" :directory '(:relative "test"))) *here*))

(defun decompose-hulls (file &rest args)
  (let* ((f (wavefront:parse file))
         (meshes (wavefront:extract-meshes f NIL '(:position)))
         (mesh (first meshes))
         (vertices (wavefront:vertex-data mesh))
         (indices (wavefront:index-data mesh)))
    (when (> (length meshes) 1)
      (break "Multiple meshes: ~A" meshes))
    (apply #'convex-covering:decompose vertices indices args)))

(defun decompose-file (in out &rest args)
  (format *trace-output* "; ~A -> ~A~%"
          (enough-namestring in *here*) (enough-namestring out *here*))
  (convex-covering::export-hulls (apply #'decompose-hulls (obj-file in) args)
                                 (merge-pathnames out (obj-file in))))

(defun decompose-test-files (&rest args)
  (let ((outdir (merge-pathnames (make-pathname :directory '(:relative "test" "output")) *here*)))
    (dolist (in (directory (obj-file (make-pathname :name "gear-1" #+no :wild))))
      (let ((out (make-pathname :name (pathname-name in) :type "obj" :defaults outdir)))
        (time (apply #'decompose-file in out args))))))

(defun test (name &key output (visualization 50))
  (let ((convex-covering::*debug-output* output)
        (convex-covering::*debug-visualizations* visualization))
    (let ((outdir (merge-pathnames (make-pathname :directory '(:relative "test" "output")) *here*)))
      (ensure-directories-exist outdir)
      (dolist (in (directory (obj-file (make-pathname :name name #+no :wild))))
        (let ((out (make-pathname :name (pathname-name in) :type "obj" :defaults outdir)))
          (time (decompose-file in out)))))))
