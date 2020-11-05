module Nav exposing (..)

import Html exposing (..)
import Html.Attributes exposing (..)


navbar =
    nav [ class "navbar navbar-expand-lg navbar-light bg-light" ]
        [ a [ class "navbar-brand text-primary", href "/" ]
            [ text "Digital Expression Explorer 2" ]
        , button
            [ attribute "aria-controls" "navbarSupportedContent"
            , attribute "aria-expanded" "false"
            , attribute "aria-label" "Toggle navigation"
            , class "navbar-toggler"
            , attribute "data-target" "#navbarSupportedContent"
            , attribute "data-toggle" "collapse"
            , type_ "button"
            ]
            [ span [ class "navbar-toggler-icon" ]
                []
            ]
        , div [ class "collapse navbar-collapse", id "navbarSupportedContent" ]
            [ ul [ class "navbar-nav mr-auto" ]
                [ li [ class "nav-item" ]
                    [ a [ class "nav-link", href "#" ]
                        [ text "Advanced Search" ]
                    ]
                , li [ class "nav-item" ]
                    [ a [ class "nav-link", href "#" ]
                        [ text "FAQ" ]
                    ]
                , li [ class "nav-item" ]
                    [ a [ class "nav-link", href "#" ]
                        [ text "About" ]
                    ]

                --, li [ class "nav-item dropdown" ]
                --    [ a [ attribute "aria-expanded" "false", attribute "aria-haspopup" "true", class "nav-link dropdown-toggle", attribute "data-toggle" "dropdown", href "#", id "navbarDropdown", attribute "role" "button" ]
                --        [ text "Dropdown        " ]
                --    , div [ attribute "aria-labelledby" "navbarDropdown", class "dropdown-menu" ]
                --        [ a [ class "dropdown-item", href "#" ]
                --            [ text "Action" ]
                --        , a [ class "dropdown-item", href "#" ]
                --            [ text "Another action" ]
                --        , div [ class "dropdown-divider" ]
                --            []
                --        , a [ class "dropdown-item", href "#" ]
                --            [ text "Something else here" ]
                --        ]
                --    ]
                ]
            ]
        ]
