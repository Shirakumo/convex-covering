(in-package #:org.shirakumo.fraf.convex-covering)

;;; `with-thread-pool' macro

(deftype thread-pool-designator ()
  '(or (eql NIL) (integer 1) lparallel:kernel (eql T)))

(defun call-with-thread-pool (thread-pool-designator continuation)
  (flet ((call-with-new-kernel (thread-count)
           (let ((kernel (lparallel:make-kernel thread-count)))
             (unwind-protect
                  (let ((lparallel:*kernel* kernel))
                    (funcall continuation))
               (let ((lparallel:*kernel* kernel))
                 (lparallel:end-kernel))))))
    (etypecase thread-pool-designator
      ((member NIL 1)
       (let ((lparallel:*kernel* nil))
         (funcall continuation)))
      ((integer 2)
       (call-with-new-kernel thread-pool-designator))
      (lparallel:kernel
       (let ((lparallel:*kernel* thread-pool-designator))
         (funcall continuation)))
      ((eql T)
       (call-with-new-kernel (org.shirakumo.machine-state:machine-cores))))))

(defmacro with-thread-pool ((thread-pool-designator) &body body)
  `(call-with-thread-pool ,thread-pool-designator (lambda () ,@body)))

;;; `maybe-plet' macro

;;; Does not support declarations.
(defmacro maybe-plet ((&rest bindings) &body body)
  `(if lparallel:*kernel*
       (lparallel:plet ,bindings ,@body)
       ,(labels ((emit (remaining)
                   (if remaining
                       (destructuring-bind (first &rest rest) remaining
                         (destructuring-bind (name-or-names value) first
                           `(multiple-value-bind ,(if (listp name-or-names)
                                                      name-or-names
                                                      (list name-or-names))
                                ,value
                              ,(emit rest))))
                       `(progn ,@body))))
          (emit bindings))))

;;; `with-tasks' macro

(defun expand-task (ordered channel task-count body)
  (let ((nid (gensym "ID")))
    `(let (,@(when ordered `((,nid ,task-count))))
       (lparallel:submit-task
        ,channel (lambda ()
                   ,@(if ordered
                         `((cons ,nid (progn ,@body)))
                         body)))
       (incf ,task-count))))

(defun expand-do-results (ordered channel task-count variable body)
  (if ordered
      (let ((nresults (gensym "RESULTS")))
        `(let* ((,nresults (loop repeat ,task-count
                                 collect (lparallel:receive-result ,channel)))
                (,nresults (sort ,nresults #'< :key #'car)))
           (mapc (lambda (,variable)
                   (let ((,variable (cdr ,variable)))
                     ,@body))
                 ,nresults)))
      `(loop repeat ,task-count
             for ,variable = (lparallel:receive-result ,channel)
             do (progn ,@body))))

(defmacro with-tasks ((&key ordered) &body body)
  (let ((nchannel (gensym "CHANNEL"))
        (ntask-count (gensym "TASK-COUNT"))
        (nresults (gensym "RESULTS")))
    `(if lparallel:*kernel*
         (let ((,nchannel (lparallel:make-channel))
               (,ntask-count 0))
           (macrolet ((task (&body body)
                        (expand-task ',ordered ',nchannel ',ntask-count body))
                      (do-results ((result) &body body)
                        (expand-do-results ',ordered ',nchannel ',ntask-count result body)))
             ,@body))
         (let ((,nresults '()))
           (macrolet ((task (&body body)
                        `(push (progn ,@body) ,',nresults))
                      (do-results ((result) &body body)
                        `(dolist (,result (nreverse ,',nresults))
                           ,@body)))
             ,@body)))))
