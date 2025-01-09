(in-package #:org.shirakumo.fraf.convex-covering)

;;; Variables

(defvar *debug-output* NIL)

(defvar *debug-visualizations* '())

(defun debug-visualizations-p (i)
  (etypecase *debug-visualizations*
    (null
     NIL)
    (integer
     (zerop (mod i *debug-visualizations*)))
    (function
     (funcall *debug-visualizations* i))))

;;; Debug output

(defun d (format-control &rest format-arguments)
  (case *debug-output*
    ((NIL))
    ((T) (apply #'format *trace-output* format-control format-arguments))
    (T (apply #'format *debug-output* format-control format-arguments))))
