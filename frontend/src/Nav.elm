module Nav exposing (..)

import HomePage.Main exposing (Msg(..))
import Html exposing (..)
import Html.Attributes exposing (..)
import Html.Events exposing (onClick)
import Themes exposing (defaultNavTheme)


navBar : Html Msg
navBar =
    let
        background =
            defaultNavTheme.background

        textColour =
            defaultNavTheme.textColour
    in
    nav ([ class "navbar navbar-expand-lg navbar-light" ] |> background |> textColour)
        [ a
            (textColour [ class "navbar-brand", href "/" ])
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
                        [ text "FAQ" ]
                    ]
                , li [ class "nav-item" ]
                    [ a [ class "nav-link", href "#" ]
                        [ text "About" ]
                    ]
                , li [ class "nav-item dropdown dropright" ]
                    [ a
                        [ attribute "aria-expanded" "false"
                        , attribute "aria-haspopup" "true"
                        , class "nav-link dropdown-toggle"
                        , attribute "data-toggle" "dropdown"
                        , href "#"
                        , id "navbarDropdown"
                        , attribute "role" "button"
                        ]
                        [ text "Search" ]
                    , div [ attribute "aria-labelledby" "navbarDropdown", class "dropdown-menu" ]
                        [ li [ class "dropdown-item", onClick SearchProjects ]
                            [ b [] [ text "Projects" ] ]
                        , li [ class "dropdown-item", onClick SearchRuns ]
                            [ b [] [ text "Runs" ] ]
                        ]
                    ]
                ]
            ]
        ]
