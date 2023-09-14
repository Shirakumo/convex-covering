(in-package #:org.shirakumo.fraf.convex-covering)

;;; Variables

(defvar *debug-output* nil)

(defvar *debug-visualizations* nil)

(defun debug-visualizations-p (i)
  (etypecase *debug-visualizations*
    (null
     nil)
    (integer
     (zerop (mod i *debug-visualizations*)))
    (function
     (funcall *debug-visualizations* i))))

;;; Debug output

(defun d (format-control &rest format-arguments)
  (when *debug-output*
    (apply #'format *trace-output* format-control format-arguments)))
